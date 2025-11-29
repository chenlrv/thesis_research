def train_dgi(data: Data, hidden_dim: int = 128, epochs: int = 200, lr: float = 1e-3,
              device: str = "cuda" if torch.cuda.is_available() else "cpu"):
    """
    Train DGI on the given graph to learn unsupervised embeddings.
    Returns the trained encoder.
    """
    x, edge_index = data.x.to(device), data.edge_index.to(device)

    encoder = GCNEncoder(in_dim=x.size(1), hidden_dim=hidden_dim)
    model = DGI(
        in_channels=x.size(1),
        hidden_channels=hidden_dim,
        encoder=encoder
    ).to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    model.train()
    for epoch in range(epochs):
        optimizer.zero_grad()
        pos_z, neg_z, summary = model(x, edge_index)
        loss = model.loss(pos_z, neg_z, summary)
        loss.backward()
        optimizer.step()

        if (epoch + 1) % 20 == 0:
            print(f"Epoch {epoch+1}/{epochs} - DGI loss: {loss.item():.4f}")

    return model.encoder.to(device)
