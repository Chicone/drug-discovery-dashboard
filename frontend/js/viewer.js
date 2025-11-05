let RDKitModuleReady = null;

window.addEventListener("DOMContentLoaded", async () => {
  console.log("Loading RDKit.js...");
  RDKitModuleReady = await initRDKitModule();
  console.log("✅ RDKit.js ready");
});

async function drawMolecule() {
  const smiles = document.getElementById("smiles").value.trim();
  if (!smiles) return alert("Enter a SMILES!");

  const container = document.getElementById("viewer3d");

  // Clear previous content (important when re-rendering)
  container.innerHTML = "";

  // Make sure container is positioned correctly
  container.style.position = "relative";

  // Fetch 3D data
  const res = await fetch(`/generate-3d?smiles=${encodeURIComponent(smiles)}`);
  const pdb = await res.text();

  if (pdb.startsWith("ERROR")) {
    alert("Invalid SMILES!");
    return;
  }

  // Create viewer *attached to the DOM element*
  const viewer = $3Dmol.createViewer(container, {
    defaultcolors: $3Dmol.rasmolElementColors,
    backgroundColor: "white",
  });

  viewer.addModel(pdb, "pdb");
  viewer.setStyle({}, { stick: { radius: 0.2, colorscheme: "Jmol" } });
  viewer.zoomTo();
  viewer.render();
}

async function loadFromFile() {
  const fileInput = document.getElementById("fileInput");
  if (!fileInput || fileInput.files.length === 0) {
    alert("Please select a file first.");
    return;
  }

  const file = fileInput.files[0];
  console.log("📄 Selected file:", file.name);

  try {
    const text = await file.text();
    console.log("📜 File content preview:", text.slice(0, 100));

    // Split into lines, ignore blanks
    const lines = text.split(/\r?\n/).filter(line => line.trim().length > 0);
    if (lines.length === 0) {
      alert("File is empty!");
      return;
    }

    // Extract first SMILES token (before space/comma/tab)
    const firstLine = lines[0].trim();
    const smiles = firstLine.split(/[ ,;\t]/)[0];
    console.log("🧪 Parsed SMILES:", smiles);

    if (!smiles || smiles.length < 2) {
      alert("No SMILES found in file.");
      return;
    }

    // Fill text box for user visibility
    document.getElementById("smiles").value = smiles;

    // Render molecule
    await drawMolecule();

  } catch (err) {
    console.error("Error reading file:", err);
    alert("Error reading file. Check console for details.");
  }
}