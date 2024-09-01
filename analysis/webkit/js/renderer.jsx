import React, { useRef, useEffect, useState } from "react";

const atom_to_color_map = {
  C: "black",
  O: "red",
  N: "blue",
  H: "white",
  Cl: "green",
  Br: "orange",
  F: "yellow",
  I: "purple",
  P: "pink",
  S: "brown",
  Si: "cyan",
  B: "magenta",
  Se: "grey",
  As: "lightblue",
  Te: "lightgreen",
  "*": "tan",
};

const Renderer = ({
  canvasRef,
  drawerEdges,
  drawerNodes,
  drawerArrows,
  idx,
}) => {
  function drawCanvas() {
    const canvas = canvasRef.current;
    const ctx = canvas.getContext("2d");
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    ctx.fillStyle = "darkgrey"; // Replace '#color' with your desired color code
    ctx.fillRect(0, 0, canvas.width, canvas.height); // if (canvasRef.current && canvasRef.current.parentNode) {

    if (
      drawerNodes === undefined ||
      drawerEdges === undefined ||
      drawerArrows === undefined
    ) {
      return;
    }
    // draw all edges

    drawerArrows.forEach((arrow) => {
      ctx.strokeStyle = "white";

      ctx.beginPath();
      ctx.moveTo(arrow.x, arrow.y);
      ctx.lineTo(arrow.x + 100, arrow.y);
      ctx.stroke();
      ctx.font = "12px Arial";

      // draw the arrow head
      ctx.beginPath();
      ctx.moveTo(arrow.x + 100, arrow.y);
      ctx.lineTo(arrow.x + 95, arrow.y + 5);
      ctx.lineTo(arrow.x + 95, arrow.y - 5);
      ctx.fillStyle = "white";
      ctx.fill();
      ctx.stroke();
    });

    drawerEdges.forEach((edge) => {
      // console.log(edge)
      const sourceNode = drawerNodes.find((node) => node.id === edge.source);
      const targetNode = drawerNodes.find((node) => node.id === edge.target);

      if (sourceNode === undefined || targetNode === undefined) {
        return;
      }

      ctx.strokeStyle = "white";

      // if edge.bondSmarts == "=" draw a double bond
      if (edge.bondSmarts === "=") {
        // console.log("HELLO:??", edge)
        ctx.beginPath();
        ctx.strokeStyle = "green";
        ctx.moveTo(sourceNode.x - 2.5, sourceNode.y - 2.5);
        ctx.lineTo(targetNode.x - 2.5, targetNode.y - 2.5);
        ctx.stroke();
        ctx.beginPath();
        ctx.moveTo(sourceNode.x + 2.5, sourceNode.y + 2.5);
        ctx.lineTo(targetNode.x + 2.5, targetNode.y + 2.5);
        ctx.stroke();
      } else if (edge.bondSmarts === "#") {
        ctx.beginPath();
        ctx.strokeStyle = "gold";
        ctx.moveTo(sourceNode.x - 3, sourceNode.y - 3);
        ctx.lineTo(targetNode.x - 3, targetNode.y - 3);
        ctx.stroke();
        ctx.beginPath();
        ctx.moveTo(sourceNode.x, sourceNode.y);
        ctx.lineTo(targetNode.x, targetNode.y);
        ctx.stroke();
        ctx.beginPath();
        ctx.moveTo(sourceNode.x + 3, sourceNode.y + 3);
        ctx.lineTo(targetNode.x + 3, targetNode.y + 3);
        ctx.stroke();
      } else {
        ctx.beginPath();
        ctx.moveTo(sourceNode.x, sourceNode.y);
        ctx.lineTo(targetNode.x, targetNode.y);
        ctx.stroke();
        // ctx.font = "12px Arial";
        // ctx.fillText(
        //   edge.bondSmarts,
        //   (sourceNode.x + targetNode.x) / 2,
        //   (sourceNode.y + targetNode.y) / 2
        // );
      }
    });

    let node_radius = 10;

    // draw all nodes
    drawerNodes.forEach((node) => {
      // console.log(node);

      ctx.beginPath();
      if (node.charge === 1) {
        ctx.strokeStyle = "red";
        ctx.lineWidth = 2;
      } else if (node.charge === -1) {
        ctx.strokeStyle = "blue";
        // change stroke width
        ctx.lineWidth = 2;
      } else {
        ctx.strokeStyle = "white";
      }
      // console.log(node.atom, node.label, node.smarts)
      if (node.smarts.includes("R")) {
        ctx.fillStyle = "white";
      } else {
        ctx.fillStyle = atom_to_color_map[node.atom];
      }
      ctx.arc(node.x, node.y, node_radius, 0, 2 * Math.PI);

      ctx.fill();
      ctx.stroke();

      const numberOfHydrogens = node.hydrogens;
      const spacing = (2 * Math.PI) / numberOfHydrogens;
      const rectWidth = 8;
      const rectHeight = 2;
      const radius = 5;

      for (let i = 0; i < numberOfHydrogens; i++) {
        const angle = i * spacing;
        const rectX =
          node.x + (radius + rectHeight / 2) * Math.cos(angle) - rectWidth / 2;
        const rectY =
          node.y + (radius + rectHeight / 2) * Math.sin(angle) - rectHeight / 2;

        ctx.save();
        ctx.translate(rectX + rectWidth / 2, rectY + rectHeight / 2);
        ctx.rotate(angle);
        ctx.fillStyle = "white";
        ctx.fillRect(-rectWidth / 2, -rectHeight / 2, rectWidth, rectHeight);
        ctx.restore();
      }

      if (node.smarts.includes("R")) {
        ctx.fillStyle = "black";
      } else {
        ctx.fillStyle = "white";
      }

      ctx.font = "12px Arial";
      let adj = 3.5;
      if (node.label.toString().length == 2) {
        adj = 7;
      }
      ctx.fillText(node.label, node.x - adj, node.y + 3.5);
    });
  }

  useEffect(() => {
    // draw all nodes and edges
    drawCanvas();
  }, [drawerNodes, drawerEdges, drawerArrows, idx]);

  useEffect(() => {
    const setCanvasSize = () => {
      if (!canvasRef.current) return; // Ensure ref is attached

      canvasRef.current.width = 750;
        // canvasRef.current.parentNode.clientWidth / 1.4285714286;
      canvasRef.current.height = 250;
    };

    // Set size initially and on window resize
    setCanvasSize();
    window.addEventListener("resize", setCanvasSize);

    // Immediate redraw once size is set
    drawCanvas(); // Ensure canvas is drawn to after resizing

    // Cleanup listener on component unmount
    return () => {
      window.removeEventListener("resize", setCanvasSize);
    };
  }, []); // Only set up and clean up once

  return <canvas ref={canvasRef} />;
};

export default Renderer;
