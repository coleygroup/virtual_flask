import React, { useRef, useEffect, useState } from "react";

function ClickableCanvas({ data, setSelected, selected }) {
  const canvasRef = useRef(null);
  const [updatedData, setUpdatedData] = useState(data);

  // Function to draw points on the canvas
  const drawPoints = (ctx, processedData) => {
    ctx.clearRect(0, 0, ctx.canvas.width, ctx.canvas.height); // Clear the canvas first
    Object.values(processedData).forEach((item) => {
      if (selected && selected.smiles === item.smiles) {
        ctx.fillStyle = "red"; // Color of the selected point
      } else {
        ctx.fillStyle = "blue"; // Color of the points
      }
      ctx.beginPath();
      ctx.arc(item.pca_x, item.pca_y, 5, 0, 2 * Math.PI); // Draw a circle for each point
      ctx.fill();
    });
  };

  // Preprocess data to fit the canvas with 10% buffer
  const preprocessData = (ctx, data) => {
    const values = Object.values(data);
    const minX = Math.min(...values.map((item) => item.pca_x));
    const maxX = Math.max(...values.map((item) => item.pca_x));
    const minY = Math.min(...values.map((item) => item.pca_y));
    const maxY = Math.max(...values.map((item) => item.pca_y));

    // Calculate ranges and factors considering a 10% padding
    const paddingX = (maxX - minX) * 0.1;
    const paddingY = (maxY - minY) * 0.1;
    const rangeX = maxX - minX + 2 * paddingX;
    const rangeY = maxY - minY + 2 * paddingY;
    const factorX = (ctx.canvas.width - 2 * paddingX) / rangeX;
    const factorY = (ctx.canvas.height - 2 * paddingY) / rangeY;
    // console.log(values);
    return values.map((item) => ({
      pca_x: (item.pca_x - minX + paddingX) * factorX,
      pca_y: (item.pca_y - minY + paddingY) * factorY,
      id: item.id,
      // data: item.data,
      // smiles: item.smiles,
      smiles: item.unmapped_smiles,
      // mapped_input_smiles: item.mapped_input_smiles,
      img: item.img,
      // mech_path: item.mech_path,
      // path_data: item.path_data,
    }));
  };

  // Handle the canvas click event
  const handleCanvasClick = (event) => {
    const rect = canvasRef.current.getBoundingClientRect();
    const x = event.clientX - rect.left; // x position within the canvas
    const y = event.clientY - rect.top; // y position within the canvas

    let closest = null;
    let minDistance = Infinity; // Start with the largest possible distance
    Object.keys(updatedData).forEach((key) => {
      const item = updatedData[key];
      const distance = Math.sqrt((x - item.pca_x) ** 2 + (y - item.pca_y) ** 2);

      if (distance < minDistance) {
        minDistance = distance; // Update the minimum distance
        closest = { key, ...item }; // Update the closest item
      }
    });

    if (closest && minDistance < 5) {
      setSelected(closest);
    } else {
      console.log("No point is within 5 pixels of the click");
    }
  };

  useEffect(() => {
    const setCanvasSize = () => {
      if (canvasRef.current && canvasRef.current.parentNode) {
        canvasRef.current.width = canvasRef.current.parentNode.clientWidth;
        canvasRef.current.height = canvasRef.current.parentNode.clientHeight;
      }
    };

    // Set size initially and on window resize
    setCanvasSize();
    window.addEventListener("resize", setCanvasSize);

    // Cleanup listener on component unmount
    return () => {
      window.removeEventListener("resize", setCanvasSize);
    };
  }, []);

  // UseEffect to draw points when data changes
  useEffect(() => {
    const canvas = canvasRef.current;
    const context = canvas.getContext("2d");

    const processedData = preprocessData(context, data);
    setUpdatedData(processedData);
    drawPoints(context, processedData);
  }, [data, selected]);

  return (
    <canvas
      ref={canvasRef}
      width={"100%"}
      height={"100%"}
      onClick={handleCanvasClick}
    />
  );
}

export default ClickableCanvas;
