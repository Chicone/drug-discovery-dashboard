from pathlib import Path

import torch
import torch.nn.functional as F
from torch_geometric.datasets import MoleculeNet
from torch_geometric.loader import DataLoader
from torch_geometric.nn import GCNConv, global_mean_pool
import torch.nn as nn


MODEL_PATH = Path("backend/services/ai/checkpoints/esol_gnn.pt")


class ESOLGNN(nn.Module):
    def __init__(self, in_channels, hidden_channels=64):
        super().__init__()
        self.conv1 = GCNConv(in_channels, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, hidden_channels)
        self.head = nn.Sequential(
            nn.Linear(hidden_channels, 32),
            nn.ReLU(),
            nn.Linear(32, 1),
        )

    def forward(self, data):
        x = self.conv1(data.x.float(), data.edge_index).relu()
        x = self.conv2(x, data.edge_index).relu()
        x = global_mean_pool(x, data.batch)
        return self.head(x).squeeze(-1)


def main():
    dataset = MoleculeNet(
        root="backend/data/ai_datasets",
        name="ESOL",
    )

    dataset = dataset.shuffle()

    n_train = int(0.8 * len(dataset))
    n_val = int(0.1 * len(dataset))

    train_dataset = dataset[:n_train]
    val_dataset = dataset[n_train:n_train + n_val]
    test_dataset = dataset[n_train + n_val:]

    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=64)
    test_loader = DataLoader(test_dataset, batch_size=64)

    in_channels = dataset.num_node_features
    model = ESOLGNN(in_channels=in_channels)

    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

    for epoch in range(1, 101):
        model.train()
        train_loss = 0.0

        for batch in train_loader:
            optimizer.zero_grad()
            pred = model(batch)
            loss = F.mse_loss(pred, batch.y.view(-1).float())
            loss.backward()
            optimizer.step()
            train_loss += loss.item() * batch.num_graphs

        train_loss /= len(train_dataset)

        if epoch % 10 == 0:
            val_rmse = evaluate(model, val_loader)
            print(
                f"Epoch {epoch:03d} | "
                f"train MSE {train_loss:.4f} | "
                f"val RMSE {val_rmse:.4f}"
            )

    test_rmse = evaluate(model, test_loader)
    print(f"Test RMSE: {test_rmse:.4f}")

    MODEL_PATH.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "model_state_dict": model.state_dict(),
            "in_channels": in_channels,
            "target": "ESOL log solubility",
        },
        MODEL_PATH,
    )

    print(f"Saved model to {MODEL_PATH}")


@torch.no_grad()
def evaluate(model, loader):
    model.eval()
    total_sq_error = 0.0
    total = 0

    for batch in loader:
        pred = model(batch)
        y = batch.y.view(-1).float()
        total_sq_error += ((pred - y) ** 2).sum().item()
        total += batch.num_graphs

    return (total_sq_error / total) ** 0.5


if __name__ == "__main__":
    main()