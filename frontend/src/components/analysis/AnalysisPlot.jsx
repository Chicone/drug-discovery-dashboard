// frontend/src/components/analysis/AnalysisPlot.jsx
import Plot from "react-plotly.js";
import { Box, Typography } from "@mui/material";


export default function AnalysisPlot({
  title = "Plot",
  plotData,
  plotMode,
  whiteBG,
  metric
}) {

  console.log("Metric:", metric);
  console.log("PlotData:", plotData);

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

      const isCom = title.includes("COM");

      return [
        {
          x: x_ns,
          y: isCom
            ? plotData.aggregate.mean.map(v => v * 10.0)
            : plotData.aggregate.mean,
          type: "scatter",
          mode: "lines",
          name: "Mean",
        },
        {
          x: x_ns,
          y: isCom
            ? plotData.aggregate.min.map(v => v * 10.0)
            : plotData.aggregate.min,
          type: "scatter",
          mode: "lines",
          name: "Min",
        },
        {
          x: x_ns,
          y: isCom
            ? plotData.aggregate.max.map(v => v * 10.0)
            : plotData.aggregate.max,
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

  // ----------------------------------------
  // ORIENTATION MODE
  // ----------------------------------------
  if (metric === "orientation" && plotData.replicas) {
  return plotData.replicas.map((rep, idx) => ({
    x: rep.times_ns,
    y: rep.angle_deg,
    type: "scatter",
    mode: "lines",
    name: rep.job_id || `Replica ${idx + 1}`,
  }));
}

  // ----------------------------------------
  // RMSD STRUCTURED RESPONSE
  // ----------------------------------------
  if (plotData.replicas) {
    return plotData.replicas.map((rep, idx) => {
      const dt = rep.timestep_ps ?? dt_ps;
      const isCom = title.includes("COM");
      const series = isCom ? rep.values : rep.rmsd;

      const x = rep.times_ps
        ? rep.times_ps.map((t) => t / 1000.0)
        : series.map((_, i) =>
            Number.isFinite(dt) ? (i * dt) / 1000.0 : i
          );

      return {
        x,
        y: isCom
          ? rep.values.map(v => v * 10.0)
          : rep.rmsd,
        type: "scatter",
        mode: "lines",
        name: rep.job_id || `Replica ${idx + 1}`,
      };
    });
  }

  // ----------------------------------------
  // COM RAW DICTIONARY
  // ----------------------------------------
  return Object.entries(plotData)
    .filter(([_, values]) => Array.isArray(values))
    .map(([jobId, values]) => ({
      x: values.map((_, i) => i),
      y: values.map((v) => v * 10.0),
      type: "scatter",
      mode: "lines",
      name: jobId,
    }));
};


const traces =
  plotMode === "aggregate"
    ? buildAggregate()
    : buildIndividual();


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
                ? `${title} — Aggregate`
                : `${title} — Individual Traces`,


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
            title: {
              text:
              metric === "orientation"
                ? "Principal Axis Angle (deg)"
                : title.includes("COM")
                ? "COM Distance (Å)"
                : "RMSD (Å)",
            },
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
