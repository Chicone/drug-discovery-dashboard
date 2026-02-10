import { Box, Typography, TextField } from "@mui/material";

export default function DistanceOptions() {
  return (
    <Box sx={{ mb: 3 }}>
      <Typography sx={{ fontWeight: 600, mb: 1 }}>Distance Measurement</Typography>

      <TextField
        label="Atom 1"
        variant="outlined"
        size="small"
        sx={{ width: "100%", mb: 1 }}
      />

      <TextField
        label="Atom 2"
        variant="outlined"
        size="small"
        sx={{ width: "100%" }}
      />
    </Box>
  );
}
