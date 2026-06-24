from pathlib import Path

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.datasets import MoleculeNet
from torch_geometric.loader import DataLoader
from torch_geometric.nn import GINConv, global_mean_pool


MODEL_PATH = Path("backend/services/ai/checkpoints/tox21_gnn.pt")

TOX21_TASKS = [
    "NR-AR",
    "NR-AR-LBD",
    "NR-AhR",
    "NR-Aromatase",
    "NR-ER",
    "NR-ER-LBD",
    "NR-PPAR-gamma",
    "SR-ARE",
    "SR-ATAD5",
    "SR-HSE",
    "SR-MMP",
    "SR-p53",
]


class Tox21GNN(nn.Module):
    def __init__(self, in_channels: int, hidden_channels: int = 64,
                 out_channels: int = 12):
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
            nn.Linear(64, out_channels),
        )

    def forward(self, data):
        x = self.conv1(data.x.float(), data.edge_index).relu()
        x = self.conv2(x, data.edge_index).relu()
        x = global_mean_pool(x, data.batch)
        return self.head(x)


def masked_bce_loss(logits, y):
    # Tox21 has missing labels marked as NaN.
    y = y.float()
    mask = ~torch.isnan(y)

    if mask.sum() == 0:
        return torch.tensor(0.0, device=logits.device)

    return F.binary_cross_entropy_with_logits(
        logits[mask],
        y[mask],
    )


@torch.no_grad()
def evaluate(model, loader):
    model.eval()
    total_loss = 0.0
    total_graphs = 0

    for batch in loader:
        logits = model(batch)
        loss = masked_bce_loss(logits, batch.y)

        total_loss += loss.item() * batch.num_graphs
        total_graphs += batch.num_graphs

    return total_loss / total_graphs


def main():
    dataset = MoleculeNet(
        root="backend/data/ai_datasets",
        name="Tox21",
    )

    dataset = dataset.shuffle()

    # define number of outputs (12 for Tox21)
    out_channels = dataset[0].y.numel()

    n_train = int(0.8 * len(dataset))
    n_val = int(0.1 * len(dataset))

    train_dataset = dataset[:n_train]
    val_dataset = dataset[n_train:n_train + n_val]
    test_dataset = dataset[n_train + n_val:]

    train_loader = DataLoader(train_dataset, batch_size=64, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=128)
    test_loader = DataLoader(test_dataset, batch_size=128)

    model = Tox21GNN(
        in_channels=dataset.num_node_features,
        out_channels=out_channels,
    )

    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

    for epoch in range(1, 51):
        model.train()
        total_loss = 0.0
        total_graphs = 0

        for batch in train_loader:
            optimizer.zero_grad()

            logits = model(batch)
            loss = masked_bce_loss(logits, batch.y)

            loss.backward()
            optimizer.step()

            total_loss += loss.item() * batch.num_graphs
            total_graphs += batch.num_graphs

        train_loss = total_loss / total_graphs

        if epoch % 5 == 0:
            val_loss = evaluate(model, val_loader)
            print(
                f"Epoch {epoch:03d} | "
                f"train BCE {train_loss:.4f} | "
                f"val BCE {val_loss:.4f}"
            )

    test_loss = evaluate(model, test_loader)
    print(f"Test BCE: {test_loss:.4f}")

    MODEL_PATH.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "model_state_dict": model.state_dict(),
            "in_channels": dataset.num_node_features,
            "out_channels": out_channels,
            "tasks": TOX21_TASKS,
            "target": "Tox21 toxicity probabilities",
        },
        MODEL_PATH,
    )

    print(f"Saved model to {MODEL_PATH}")


if __name__ == "__main__":
    main()