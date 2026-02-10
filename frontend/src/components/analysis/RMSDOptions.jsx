import { Box, Typography, TextField } from "@mui/material";

export default function RMSDOptions() {
  return (
    <Box sx={{ mb: 3 }}>
      <Typography sx={{ fontWeight: 600, mb: 1 }}>RMSD Settings</Typography>

      <TextField
        label="Reference Frame"
        type="number"
        variant="outlined"
        size="small"
        sx={{
          mt: 1,
          width: "100%",
          input: { color: "white" },
          label: { color: "#aaa" },
        }}
      />
    </Box>
  );
}
