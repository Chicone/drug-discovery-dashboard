import { useEffect, useRef } from "react";

function MoleculeViewer({ smiles }) {
  const viewerRef = useRef(null);

  useEffect(() => {
    if (!smiles || !window.$3Dmol) return;

    const element = viewerRef.current;
    element.innerHTML = "";

    // Explicitly set positioning context before creating viewer
    element.style.position = "relative";

    const viewer = window.$3Dmol.createViewer(element, {
      backgroundColor: "#111111",
    });

    fetch(`/api/mol3d?smiles=${encodeURIComponent(smiles)}`)
      .then((r) => r.text())
      .then((molData) => {
        viewer.addModel(molData, "mol");
        viewer.setStyle({}, { stick: {}, sphere: { scale: 0.3 } });
        viewer.zoomTo();
        viewer.render();
      });
  }, [smiles]);

  return (
    <div
      ref={viewerRef}
      style={{
        width: "600px",
        height: "600px",
        border: "1px solid #555",
        borderRadius: "6px",
        background: "#111",
        marginTop: "1rem",
        position: "relative", // <---- THIS IS THE FIX
        overflow: "hidden",   // ensure canvas stays inside
      }}
    />
  );
}

export default MoleculeViewer;

