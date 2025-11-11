import { motion, AnimatePresence } from "framer-motion";
import { useState } from "react";
import {
  AppBar,
  Toolbar,
  Typography,
  Tabs,
  Tab,
  Box,
  Container,
  GlobalStyles,
} from "@mui/material";
import { useTheme } from "@mui/material/styles";

import MolecularDesign from "./pages/MolecularDesign.jsx";
import Docking from "./pages/Docking.jsx";
import AdmetAnalysis from "./pages/AdmetAnalysis.jsx";
import Library from "./pages/Library.jsx";
import Reports from "./pages/Reports.jsx";
import MolecularDynamics from "./pages/MolecularDynamics.jsx";

import logo from "./assets/logo.png";

export default function App() {
  const theme = useTheme();
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
    <Box
      sx={{
        background: "#121212",
        color: "white",
        minHeight: "100vh",
        overflowX: "hidden",
      }}
    >
      <GlobalStyles
        styles={{
          html: { scrollbarGutter: "stable both-edges" },
          body: {
            scrollbarGutter: "stable both-edges",
            overflowX: "hidden !important",
          },
          "#root": {
            overflowX: "hidden !important",
          },
        }}
      />

      <AppBar
        position="sticky"
        sx={{
          background: "rgba(18, 18, 18, 0.8)",
          backdropFilter: "blur(8px)",
          borderBottom: "1px solid rgba(255,255,255,0.1)",
        }}
      >
        <Toolbar>
          {/* LEFT SIDE: LOGO + TITLE */}
          <Box
            sx={{
              display: "flex",
              alignItems: "center",
              flexGrow: 1,
              minWidth: 0,
            }}
          >
            <Box
              component="img"
              src={logo}
              alt="Logo"
              sx={{
                width: 36,
                height: 36,
                mr: 1.5,
                objectFit: "contain",
                filter: "drop-shadow(0 0 2px rgba(0,0,0,0.3))",
              }}
            />
            <Box>
              <Typography
                variant="h6"
                sx={{
                  fontWeight: 600,
                  color: theme.palette.primary.main,
                  lineHeight: 1.2,
                }}
              >
                Drug Design Dashboard
              </Typography>
              <Typography
                variant="caption"
                sx={{ color: "text.secondary", letterSpacing: 1.5 }}
              >
                UNIGE · Pharmaceutical Sciences
              </Typography>
            </Box>
          </Box>
          <Box sx={{ flexGrow: 1 }} />
          <Box sx={{ width: 40 }} />
          {/* RIGHT SIDE: NAVIGATION TABS */}
          <Tabs
            value={tab}
            onChange={(_, v) => setTab(v)}
            textColor="inherit"
            indicatorColor="secondary"
            variant="scrollable"
            scrollButtons="auto"
            sx={{
              flexShrink: 0,
              "& .MuiTabs-flexContainer": {
                justifyContent: "flex-end",
              },
              "& .MuiTab-root": {
                minHeight: 64,
                textTransform: "none",
                fontWeight: 500,
                whiteSpace: "nowrap",
              },
            }}
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

      {/* ✅ CONTENT SECTION */}
      <Container
        maxWidth="lg"
        sx={{
          py: 4,
          minHeight: "calc(100vh - 64px)",
          display: "flex",
          flexDirection: "column",
          overflowX: "hidden",
        }}
      >
        <AnimatePresence mode="wait">
          <motion.div
            key={tab}
            initial={{ opacity: 0, y: 8 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: -8 }}
            transition={{ duration: 0.3 }}
          >
            {render()}
          </motion.div>
        </AnimatePresence>
      </Container>
    </Box>
  );
}
