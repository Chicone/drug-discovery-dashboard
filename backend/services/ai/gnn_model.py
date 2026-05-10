import torch
import torch.nn as nn
from torch_geometric.nn import GCNConv, global_mean_pool


class SimpleMolecularGNN(nn.Module):
    def __init__(self, in_channels: int = 4, hidden_channels: int = 64):
        super().__init__()

        # Graph convolution layers (message passing)
        self.conv1 = GCNConv(in_channels, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, hidden_channels)

        # Final prediction head
        self.head = nn.Sequential(
            nn.Linear(hidden_channels, 32),
            nn.ReLU(),
            nn.Linear(32, 1),
        )

    def forward(self, data):
        x = data.x
        edge_index = data.edge_index

        # Propagate information across graph
        x = self.conv1(x, edge_index).relu()
        x = self.conv2(x, edge_index).relu()

        # Create batch (single molecule → all zeros)
        if hasattr(data, "batch") and data.batch is not None:
            batch = data.batch
        else:
            batch = torch.zeros(x.size(0), dtype=torch.long, device=x.device)

        # Pool node features → molecule-level vector
        pooled = global_mean_pool(x, batch)

        # Predict scalar property
        return self.head(pooled)