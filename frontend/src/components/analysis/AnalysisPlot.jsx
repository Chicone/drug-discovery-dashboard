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
  // Build X-axis depending on metric
  let x_ns = [];

  if (metric === "activation") {
    // Use the length of tm3_tm6_distance.mean
    const n = plotData?.aggregate?.tm3_tm6_distance?.mean?.length ?? 0;
    x_ns = [...Array(n).keys()].map(i =>
      Number.isFinite(dt_ps) ? (i * dt_ps) / 1000.0 : i
    );
  } else {
    // RMSD and COM use generic aggregate.mean
    x_ns =
      plotData?.aggregate?.mean?.map((_, i) =>
        Number.isFinite(dt_ps) ? (i * dt_ps) / 1000.0 : i
      ) ?? [];
  }
  // -----------------------------------------------------
  // Build AGGREGATE traces
  // -----------------------------------------------------
    const buildAggregate = () => {
      if (!plotData.aggregate) return [];

      // ----------------------------------------
      // ACTIVATION METRICS (override)
      // ----------------------------------------
      if (metric === "activation") {
        const exp = v => v * 10.0;  // nm → Å

        const dist = plotData.aggregate.tm3_tm6_distance;
        const disp = plotData.aggregate.tm6_displacement;

        return [
          // TM3–TM6 distance
          {
            x: x_ns,
            y: dist.mean.map(exp),
            type: "scatter",
            mode: "lines",
            name: "TM3–TM6 mean (Å)",
          },
          {
            x: x_ns,
            y: dist.min.map(exp),
            type: "scatter",
            mode: "lines",
            name: "TM3–TM6 min (Å)",
          },
          {
            x: x_ns,
            y: dist.max.map(exp),
            type: "scatter",
            mode: "lines",
            name: "TM3–TM6 max (Å)",
          },

          // TM6 displacement along axis
          {
            x: x_ns,
            y: disp.mean.map(exp),
            type: "scatter",
            mode: "lines",
            name: "TM6 disp mean (Å)",
          },
          {
            x: x_ns,
            y: disp.min.map(exp),
            type: "scatter",
            mode: "lines",
            name: "TM6 disp min (Å)",
          },
          {
            x: x_ns,
            y: disp.max.map(exp),
            type: "scatter",
            mode: "lines",
            name: "TM6 disp max (Å)",
          },
        ];
      }

      // ----------------------------------------
      // GENERIC MODE (RMSD / COM)
      // ----------------------------------------
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
      // ACTIVATION METRICS
      // ----------------------------------------
      if (metric === "activation" && plotData.replicas) {

        // TM3–TM6 IC distance
        const distTraces = plotData.replicas.map((rep, idx) => {
          const dt = rep.timestep_ps ?? dt_ps;

          const x = rep.tm3_tm6_distance.map(
            (_, i) => (i * dt) / 1000.0
          );

          return {
            x,
            y: rep.tm3_tm6_distance.map(v => v * 10.0), // nm → Å
            type: "scatter",
            mode: "lines",
            name: (rep.job_id || rep.name || "").slice(0, 8) + "…" + " (TM3–TM6)",
          };
        });

        // TM6 displacement
        const dispTraces = plotData.replicas.map((rep, idx) => {
          const dt = rep.timestep_ps ?? dt_ps;

          const x = rep.tm6_displacement.map(
            (_, i) => (i * dt) / 1000.0
          );

          return {
            x,
            y: rep.tm6_displacement.map(v => v * 10.0),
            type: "scatter",
            mode: "lines",
            name: (rep.job_id || rep.name || "").slice(0, 8) + "…" + " (TM6 disp)",
          };
        });

        return [...distTraces, ...dispTraces];
      }

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
                : metric === "activation"
                ? "TM3–TM6 dist / TM6 disp (Å)"
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
