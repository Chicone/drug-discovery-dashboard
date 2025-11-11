import { useState } from "react";
import {
  AppBar,
  Toolbar,
  Typography,
  Tabs,
  Tab,
  Box,
  Container,
} from "@mui/material";

// ✅ Make sure every one of these files exists:
import MolecularDesign from "./pages/MolecularDesign.jsx";
import Docking from "./pages/Docking.jsx";
import AdmetAnalysis from "./pages/AdmetAnalysis.jsx";
import Library from "./pages/Library.jsx";
import Reports from "./pages/Reports.jsx";
import MolecularDynamics from "./pages/MolecularDynamics.jsx";

export default function App() {
  const [tab, setTab] = useState(0);

  const render = () => {
    switch (tab) {
      case 0:
        return <MolecularDesign />;
      case 1:
        return <Docking />;
      case 2:
        return <AdmetAnalysis />;
      case 3:
        return <Library />;
      case 4:
        return <Reports />;
      case 5:
        return <MolecularDynamics />;
      default:
        return null;
    }
  };

  return (
    <Box sx={{ background: "#121212", color: "white", minHeight: "100vh" }}>
      <AppBar position="static" sx={{ background: "#1f1f1f" }}>
        <Toolbar>
          <Typography variant="h6" sx={{ flexGrow: 1 }}>
            🧬 Drug Discovery Dashboard
          </Typography>
          <Tabs
            value={tab}
            onChange={(_, v) => setTab(v)}
            textColor="inherit"
            indicatorColor="secondary"
          >
            <Tab label="Molecular Design" />
            <Tab label="Docking" />
            <Tab label="ADMET" />
            <Tab label="Library" />
            <Tab label="Reports" />
            <Tab label="Molecular Dynamics" />
          </Tabs>
        </Toolbar>
      </AppBar>

      <Container sx={{ py: 4 }}>{render()}</Container>
    </Box>
  );
}
