import { useState, useEffect } from "react";
import MoleculeViewer from "../MoleculeViewer"; // updated relative path

function MolecularDesign() {
  const [smiles, setSmiles] = useState("");
  const [data, setData] = useState(null);
  const [error, setError] = useState(null);
  const [RDKit, setRDKit] = useState(null);

  useEffect(() => {
    const loadRDKit = async () => {
      const RDKitModule = await window.initRDKitModule();
      setRDKit(RDKitModule);
    };
    if (!window.initRDKitModule) {
      // dynamically load RDKit.js script
      const script = document.createElement("script");
      script.src = "https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js";
      script.onload = loadRDKit;
      document.body.appendChild(script);
    } else {
      loadRDKit();
    }
  }, []);

  const fetchProperties = async () => {
    setError(null);
    try {
      const res = await fetch(`/api/properties?smiles=${encodeURIComponent(smiles)}`);
      const json = await res.json();
      if (json.error) setError(json.error);
      else setData(json);
    } catch (e) {
      setError("Could not connect to backend");
    }
  };

  const renderMolecule = () => {
    if (!RDKit || !smiles) return null;
    try {
      const mol = RDKit.get_mol(smiles);
      const svg = mol.get_svg();
      mol.delete();
      return (
        <div
          dangerouslySetInnerHTML={{ __html: svg }}
          style={{ border: "1px solid #ccc", marginTop: "1rem" }}
        />
      );
    } catch {
      return <p style={{ color: "red" }}>Invalid SMILES for visualization</p>;
    }
  };

return (
  <div style={{
    padding: "2rem",
    fontFamily: "sans-serif",
    color: "white",
    backgroundColor: "#121212",
    minHeight: "100vh",
    display: "flex",
    flexDirection: "column",
    alignItems: "center",
  }}>
    <h1 style={{ marginBottom: "1rem" }}>Drug Discovery Dashboard</h1>

    <div style={{ display: "flex", gap: "0.5rem" }}>
      <input
        type="text"
        placeholder="Enter SMILES (e.g. CCO)"
        value={smiles}
        onChange={(e) => setSmiles(e.target.value)}
        style={{
          width: "300px",
          padding: "0.5rem",
          borderRadius: "4px",
          border: "1px solid #555",
          color: "#fff",
          background: "#222",
        }}
      />
      <input
          type="file"
          accept=".smi,.txt"
          onChange={(e) => {
            const file = e.target.files[0];
            if (!file) return;
            const reader = new FileReader();
            reader.onload = (ev) => {
              const text = ev.target.result;
              const firstLine = text.split("\n").find((l) => l.trim().length > 0);
              if (firstLine) {
                const smilesFromFile = firstLine.split(/\s+/)[0]; // get first token
                setSmiles(smilesFromFile);
              }
            };
            reader.readAsText(file);
          }}
          style={{ marginTop: "1rem", display: "block" }}
        />
      <button
        onClick={fetchProperties}
        style={{
          padding: "0.5rem 1rem",
          background: "#444",
          border: "none",
          color: "white",
          borderRadius: "4px",
          cursor: "pointer",
        }}
      >
        Analyze
      </button>
    </div>

    {error && <p style={{ color: "red" }}>{error}</p>}

    {data && (
      <div style={{ marginTop: "2rem", textAlign: "center" }}>
        <table
          border="1"
          cellPadding="5"
          style={{
            margin: "0 auto",
            borderCollapse: "collapse",
            minWidth: "250px",
          }}
        >
          <tbody>
            {Object.entries(data).map(([key, value]) => (
              <tr key={key}>
                <td><strong>{key}</strong></td>
                <td>{value}</td>
              </tr>
            ))}
          </tbody>
        </table>

        {/* Viewer container */}
        <div style={{ marginTop: "1.5rem", display: "flex", justifyContent: "center" }}>
          <MoleculeViewer smiles={smiles} />
        </div>
      </div>
    )}
  </div>
);
}

export default MolecularDesign;
