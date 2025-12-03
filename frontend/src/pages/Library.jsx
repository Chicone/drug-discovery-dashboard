import { Paper, Typography, Table, TableBody, TableRow, TableCell } from "@mui/material";

function Library() {
  const compounds = [
    { name: "Compound A", logP: 2.5, activity: "Active" },
    { name: "Compound B", logP: 1.8, activity: "Inactive" },
    { name: "Compound C", logP: 3.1, activity: "Active" },
  ];

  return (
    <Paper sx={{ p: 3, background: "#1e1e1e" }}>
      <Typography variant="h5" gutterBottom>
        📚 Compound Library
      </Typography>
      <Typography variant="body1" sx={{ mb: 2 }}>
        Manage and review candidate molecules with their computed or predicted properties.
      </Typography>

      <Table>
        <TableBody>
          {compounds.map((c, i) => (
            <TableRow key={i}>
              <TableCell sx={{ color: "white" }}>{c.name}</TableCell>
              <TableCell sx={{ color: "white" }}>logP: {c.logP}</TableCell>
              <TableCell sx={{ color: "white" }}>{c.activity}</TableCell>
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </Paper>
  );
}

export default Library;
