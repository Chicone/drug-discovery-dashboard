from pathlib import Path

import torch
import torch.nn as nn
import torch.nn.functional as F
from tdc.single_pred import Tox
from torch_geometric.loader import DataLoader
from torch_geometric.nn import GINConv, global_mean_pool
from torch_geometric.utils.smiles import from_smiles
from pathlib import Path

MODEL_PATH = Path("backend/services/ai/checkpoints/ld50_gnn.pt")

class LD50GNN(nn.Module):
    def __init__(self, in_channels: int, hidden_channels: int = 64):
        super().__init__()

        self.conv1 = GINConv(
            nn.Sequential(
                nn.Linear(in_channels, hidden_channels),
                nn.ReLU(),
                nn.Linear(hidden_channels, hidden_channels),
            )
        )

        self.conv2 = GINConv(
            nn.Sequential(
                nn.Linear(hidden_channels, hidden_channels),
                nn.ReLU(),
                nn.Linear(hidden_channels, hidden_channels),
            )
        )

        self.head = nn.Sequential(
            nn.Linear(hidden_channels, 64),
            nn.ReLU(),
            nn.Linear(64, 1),
        )

    def forward(self, data):
        x = self.conv1(data.x.float(), data.edge_index).relu()
        x = self.conv2(x, data.edge_index).relu()
        x = global_mean_pool(x, data.batch)
        return self.head(x).squeeze(-1)

def dataframe_to_graphs(df):
    graphs = []
    failed = 0

    for _, row in df.iterrows():
        smiles = row["Drug"]
        y = float(row["Y"])

        try:
            graph = from_smiles(smiles)

            if graph.x is None or graph.x.numel() == 0:
                failed += 1
                continue

            graph.y = torch.tensor([y], dtype=torch.float)
            graphs.append(graph)

        except Exception as exc:
            failed += 1
            if failed <= 5:
                print(f"[WARN] Failed SMILES: {smiles} | {exc}")

    print(f"Converted {len(graphs)} molecules, failed {failed}")
    return graphs


@torch.no_grad()
def evaluate(model, loader):
    model.eval()
    total_abs_error = 0.0
    total = 0

    for batch in loader:
        pred = model(batch)
        y = batch.y.view(-1).float()

        total_abs_error += torch.abs(pred - y).sum().item()
        total += batch.num_graphs

    return total_abs_error / total


def main():
    data = Tox(name="LD50_Zhu")
    split = data.get_split(method="scaffold")

    train_dataset = dataframe_to_graphs(split["train"])
    val_dataset = dataframe_to_graphs(split["valid"])
    test_dataset = dataframe_to_graphs(split["test"])

    train_loader = DataLoader(train_dataset, batch_size=64, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=128)
    test_loader = DataLoader(test_dataset, batch_size=128)

    in_channels = train_dataset[0].x.shape[1]
    model = LD50GNN(in_channels=in_channels)

    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

    for epoch in range(1, 51):
        model.train()
        total_loss = 0.0
        total_graphs = 0

        for batch in train_loader:
            optimizer.zero_grad()

            pred = model(batch)
            y = batch.y.view(-1).float()

            loss = F.mse_loss(pred, y)
            loss.backward()
            optimizer.step()

            total_loss += loss.item() * batch.num_graphs
            total_graphs += batch.num_graphs

        train_mse = total_loss / total_graphs

        if epoch % 5 == 0:
            val_mae = evaluate(model, val_loader)
            print(
                f"Epoch {epoch:03d} | "
                f"train MSE {train_mse:.4f} | "
                f"val MAE {val_mae:.4f}"
            )

    test_mae = evaluate(model, test_loader)
    print(f"Test MAE: {test_mae:.4f}")

    MODEL_PATH.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "model_state_dict": model.state_dict(),
            "in_channels": in_channels,
            "target": "LD50 acute toxicity",
            "unit": "log(1/(mol/kg))",
        },
        MODEL_PATH,
    )

    print(f"Saved model to {MODEL_PATH}")


if __name__ == "__main__":
    main()