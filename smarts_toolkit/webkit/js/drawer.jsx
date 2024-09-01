import React, { useRef, useEffect, useState } from "react";

const Drawer = ({
  canvasRef,
  setDrawerNodes,
  drawerNodes,
  setDrawerEdges,
  drawerEdges,
  drawMode,
  atomSmarts,
  setAtomSmarts,
  setSelectedAtom,
  selectedAtom,
  setSelectedBond,
  selectedBond,
  setNodeId,
  nodeId,
  selectedAtoms,
  setSelectedAtoms,
  bondSmarts,
  bondId,
  setBondId,
  currAtomMap,
  setCurrAtomMap,
}) => {
  const [isDrawing, setIsDrawing] = useState(false);
  const [lastX, setLastX] = useState(0);
  const [lastY, setLastY] = useState(0);

  useEffect(() => {
    const canvas = canvasRef.current;
    const handleClick = (event) => {
      if (drawMode == "atom") {
        const nodeRadius = 35;
        const clickedNode = drawerNodes.find((node) => {
          return (
            Math.sqrt(
              (event.offsetX - node.x) ** 2 + (event.offsetY - node.y) ** 2
            ) < nodeRadius
          );
        });

        if (clickedNode) {
          setSelectedAtom(clickedNode);
          return;
        }

        setSelectedAtom(null);

        setDrawerNodes((prev) => [
          ...prev,
          {
            x: event.offsetX,
            y: event.offsetY,
            atomSmarts: atomSmarts,
            id: nodeId,
          },
        ]);

        setCurrAtomMap(currAtomMap + 1);
        let num = currAtomMap + 1;
        setAtomSmarts("[C;H0;+0:" + num + "]");

        setNodeId((prev) => prev + 1);
      } else if (drawMode == "bond") {
        if (selectedAtoms.length == 1) {
          // check if bond exists
          const nodeRadius = 35;
          const clickedNode = drawerNodes.find((node) => {
            return (
              Math.sqrt(
                (event.offsetX - node.x) ** 2 + (event.offsetY - node.y) ** 2
              ) < nodeRadius
            );
          });

          if (!clickedNode) {
            return;
          }

          for (let edge of drawerEdges) {
            if (
              edge.source == selectedAtoms[0].id &&
              edge.target == clickedNode.id
            ) {
              setSelectedAtoms([]);
              setSelectedBond(edge);
              return;
            } else if (
              edge.source == clickedNode.id &&
              edge.target == selectedAtoms[0].id
            ) {
              setSelectedAtoms([]);
              setSelectedBond(edge);
              return;
            }
          }

          setDrawerEdges((prev) => [
            ...prev,
            {
              source: selectedAtoms[0].id,
              target: clickedNode.id,
              bondSmarts: bondSmarts,
              id: bondId,
            },
          ]);
          setBondId((prev) => prev + 1);
          setSelectedAtoms([]);
        } else {
          const nodeRadius = 35;
          const clickedNode = drawerNodes.find((node) => {
            return (
              Math.sqrt(
                (event.offsetX - node.x) ** 2 + (event.offsetY - node.y) ** 2
              ) < nodeRadius
            );
          });

          if (clickedNode) {
            setSelectedAtoms((prev) => [...prev, clickedNode]);
            return;
          }
        }
      }
    };

    // if (isDragging) {
    //   return;
    // }
    canvas.addEventListener("click", handleClick);

    return () => canvas.removeEventListener("click", handleClick);
  }, [atomSmarts, drawerNodes, nodeId, bondId, drawMode, selectedAtoms]);

  useEffect(() => {
    // draw all nodes and edges
    const canvas = canvasRef.current;
    const ctx = canvas.getContext("2d");
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    // draw all edges

    drawerEdges.forEach((edge) => {
      const sourceNode = drawerNodes.find((node) => node.id === edge.source);
      const targetNode = drawerNodes.find((node) => node.id === edge.target);

      if (drawMode == "bond") {
        if (selectedBond === null) {
          ctx.strokeStyle = "black";
        } else {
          if (edge.id === selectedBond.id) {
            ctx.strokeStyle = "red";
          } else {
            ctx.strokeStyle = "black";
          }
        }
      } else {
        ctx.strokeStyle = "black";
      }

      ctx.beginPath();
      ctx.moveTo(sourceNode.x, sourceNode.y);
      ctx.lineTo(targetNode.x, targetNode.y);
      ctx.stroke();
      ctx.font = "12px Arial";
      ctx.fillText(
        edge.bondSmarts,
        (sourceNode.x + targetNode.x) / 2,
        (sourceNode.y + targetNode.y) / 2
      );
    });

    // draw all nodes
    drawerNodes.forEach((node) => {
      ctx.beginPath();

      if (drawMode == "atom") {
        if (selectedAtom === null) {
          ctx.strokeStyle = "black";
        } else {
          if (node.id === selectedAtom.id) {
            ctx.strokeStyle = "red";
          } else {
            ctx.strokeStyle = "black";
          }
        }
      } else if (drawMode == "bond") {
        if (selectedAtoms.length == 0) {
          ctx.strokeStyle = "black";
        } else {
          if (selectedAtoms.length == 1) {
            if (node.id === selectedAtoms[0].id) {
              ctx.strokeStyle = "green";
            } else {
              ctx.strokeStyle = "black";
            }
          } else if (selectedAtoms.length == 2) {
            if (
              node.id === selectedAtoms[0].id ||
              node.id === selectedAtoms[1].id
            ) {
              ctx.strokeStyle = "green";
            } else {
              ctx.strokeStyle = "black";
            }
          } else {
            ctx.strokeStyle = "black";
          }
        }
      }
      ctx.fillStyle = "white";
      ctx.arc(node.x, node.y, 35, 0, 2 * Math.PI);
      ctx.fill();
      ctx.stroke();
      ctx.fillStyle = "black";
      ctx.font = "12px Arial";
      ctx.fillText(node.atomSmarts, node.x - 31, node.y + 4);
    });
  }, [
    selectedAtom,
    drawerNodes,
    drawMode,
    selectedAtoms,
    selectedBond,
    drawerEdges,
  ]);

  const [isDragging, setIsDragging] = useState(false);
  const [selectedNodeIndex, setSelectedNodeIndex] = useState(null);
  const handleMouseDown = (e) => {
    const rect = canvasRef.current.getBoundingClientRect();
    const mouseX = e.clientX - rect.left;
    const mouseY = e.clientY - rect.top;
    const clickedNodeIndex = drawerNodes.findIndex(
      (node) =>
        Math.sqrt((node.x - mouseX) ** 2 + (node.y - mouseY) ** 2) < 35
    );
    console.log("hello!", mouseX, mouseY, clickedNodeIndex)
    for (let i = 0; i < drawerNodes.length; i++) {
      console.log(drawerNodes[i].x, drawerNodes[i].y);
    }

    if (clickedNodeIndex !== -1) {
      setIsDragging(true);
      setSelectedNodeIndex(clickedNodeIndex);
    }
  };

  const handleMouseUp = () => {
    setIsDragging(false);
    setSelectedNodeIndex(null);
  };

  const handleMouseMove = (e) => {
    if (isDragging && selectedNodeIndex !== null) {
      const rect = canvasRef.current.getBoundingClientRect();
      const newX = e.clientX - rect.left;
      const newY = e.clientY - rect.top;
      const updatedNodes = drawerNodes.map((node, index) => {
        if (index === selectedNodeIndex) {
          return { ...node, x: newX, y: newY };
        }
        return node;
      });

      setDrawerNodes(updatedNodes);
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
    window.addEventListener('resize', setCanvasSize);

    // Cleanup listener on component unmount
    return () => {
      window.removeEventListener('resize', setCanvasSize);
    };
  }, []);

  return (
    <canvas
      ref={canvasRef}
      onMouseDown={handleMouseDown}
      onMouseUp={handleMouseUp}
      onMouseOut={handleMouseUp}
      onMouseMove={handleMouseMove}
    />
  );
};

export default Drawer;
