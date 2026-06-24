import { Box, Typography } from "@mui/material";

export default function AnalysisSidebar({ children }) {
  return (
    <Box
      sx={{
        width: 280,
        p: 2,
        borderRight: "1px solid rgba(255,255,255,0.1)",
        background: "#1a1a1a",
      }}
    >
      <Typography variant="h6" sx={{ mb: 2 }}>
        Analysis Tools
      </Typography>

      {children}
    </Box>
  );
}
