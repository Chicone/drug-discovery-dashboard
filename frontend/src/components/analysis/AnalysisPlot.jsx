// frontend/src/components/analysis/AnalysisPlot.jsx
import Plot from "react-plotly.js";
import { Box, Typography } from "@mui/material";


export default function AnalysisPlot({ title = "Plot", plotData, plotMode, whiteBG }) {
  // When no data yet
  if (!plotData) {
    return (
      <Box sx={{ flex: 1, p: 2 }}>
        <Typography variant="h6" sx={{ mb: 2 }}>
          {title}
        </Typography>
        <Box
          sx={{
            height: 500,
            border: "1px solid rgba(255,255,255,0.1)",
            borderRadius: 2,
            background: "#111",
            p: 2,
            color: "gray",
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
          }}
        >
          Select jobs and run analysis.
        </Box>
      </Box>
    );
  }

  // -----------------------------------------------------
  // Extract timestep (ps)
  // -----------------------------------------------------
  const dt_ps = Number(
    plotData?.replicas?.find((r) => r?.timestep_ps != null)?.timestep_ps
  );

  // -----------------------------------------------------
  // X-axis for AGGREGATE plot
  // -----------------------------------------------------
  const x_ns =
    plotData?.aggregate?.mean?.map((_, i) =>
      Number.isFinite(dt_ps) ? (i * dt_ps) / 1000.0 : i
    ) ?? [];

  // -----------------------------------------------------
  // Build AGGREGATE traces
  // -----------------------------------------------------
  const buildAggregate = () => {
    if (!plotData.aggregate) return [];

    return [
      {
        x: x_ns,
        y: plotData.aggregate.mean,
        type: "scatter",
        mode: "lines",
        name: "Mean RMSD",
      },
      {
        x: x_ns,
        y: plotData.aggregate.min,
        type: "scatter",
        mode: "lines",
        name: "Min",
      },
      {
        x: x_ns,
        y: plotData.aggregate.max,
        type: "scatter",
        mode: "lines",
        name: "Max",
      },
    ];
  };

  // -----------------------------------------------------
  // Build INDIVIDUAL traces (one per job)
  // -----------------------------------------------------
  const buildIndividual = () => {
    if (!plotData.replicas) return [];

    return plotData.replicas.map((rep, idx) => {
      const dt = rep.timestep_ps ?? dt_ps;
      const x = rep.times_ps
      ? rep.times_ps.map(t => t / 1000.0)
      : rep.rmsd.map((_, i) => Number.isFinite(dt) ? (i * dt) / 1000.0 : i);


      return {
        x,
        y: rep.rmsd,
        type: "scatter",
        mode: "lines",
        name: rep.job_id || `Replica ${idx + 1}`,
      };
    });
  };

  const traces =
    plotMode === "aggregate" ? buildAggregate() : buildIndividual();

  return (
    <Box sx={{ flex: 1, p: 2 }}>
      <Typography variant="h6" sx={{ mb: 2 }}>
        {title}
      </Typography>

      <Box
        sx={{
          height: 500,
          border: "1px solid rgba(255,255,255,0.1)",
          borderRadius: 2,
          background: "#111",
          p: 2,
        }}
      >
        <Plot
          data={traces}
          layout={{
          title:
            plotMode === "aggregate"
              ? "Ligand RMSD — Aggregate"
              : "Ligand RMSD — Individual Traces",

          // Switch styles depending on toggle:
          template: whiteBG ? "plotly_white" : undefined,
          paper_bgcolor: whiteBG ? "white" : "#111",
          plot_bgcolor: whiteBG ? "white" : "#111",

          font: { color: whiteBG ? "black" : "white" },

          margin: { t: 60, r: 40, b: 80, l: 60 },

          xaxis: {
            title: { text: "Time (ns)" },
            showline: true,
            linecolor: whiteBG ? "black" : "white",
          },

          yaxis: {
            title: { text: "RMSD (Å)" },
            showline: true,
            linecolor: whiteBG ? "black" : "white",
          },
        }}

          style={{ width: "100%", height: "100%" }}
          config={{ displayModeBar: true }}
        />
      </Box>
    </Box>
  );
}
