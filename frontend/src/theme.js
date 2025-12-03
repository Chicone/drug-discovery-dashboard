import { createTheme } from "@mui/material/styles";

const theme = createTheme({
  palette: {
    mode: "dark",
    background: { default: "#121212", paper: "#1e1e1e" },
    primary: { main: "#00bfa5" },
    secondary: { main: "#00e676" },
    text: { primary: "#ffffff", secondary: "#aaaaaa" },
  },
  typography: {
    fontFamily: "'Inter', 'Roboto', 'Helvetica', 'Arial', sans-serif",
    h4: {
      fontWeight: 600,
      borderBottom: "2px solid #00bfa5",
      paddingBottom: "0.3rem",
      display: "inline-block",
      marginBottom: "1rem",
    },
    h5: { fontWeight: 600 },
    body1: { lineHeight: 1.6 },
    button: { textTransform: "none", fontWeight: 500 },
  },
  components: {
    MuiPaper: {
      styleOverrides: {
        root: {
          borderRadius: 10,
          padding: "1.5rem",
          transition: "0.2s",
          "&:hover": { boxShadow: "0 4px 12px rgba(0,0,0,0.3)" },
        },
      },
    },
    MuiAppBar: {
      styleOverrides: { root: { boxShadow: "none", borderBottom: "1px solid #333" } },
    },
    MuiTabs: {
      styleOverrides: { indicator: { height: 3, borderRadius: 3 } },
    },
    MuiTab: {
      styleOverrides: {
        root: {
          textTransform: "none",
          fontWeight: 500,
          fontSize: "0.95rem",
          "&:hover": { color: "#00e676" },
          "&.Mui-selected": { color: "#00bfa5" },
        },
      },
    },
    MuiButton: {
      styleOverrides: {
        root: {
          borderRadius: 8,
          fontWeight: 500,
          padding: "0.5rem 1rem",
        },
      },
    },
  },
});

export default theme;
