import React, { useState, useEffect, useRef } from "react";
import Drawer from "./drawer";

const DrawerView = () => {
  const canvasDrawerRef = useRef(null);

  const [drawerNodes, setDrawerNodes] = useState([]);
  const [drawerEdges, setDrawerEdges] = useState([]);
  const [drawMode, setDrawMode] = useState("atom");
  const [atomSmarts, setAtomSmarts] = useState("[C;H2;+0:1]");
  const [currAtomMap, setCurrAtomMap] = useState(1);
  const [selectedAtom, setSelectedAtom] = useState(null);
  const [selectedBond, setSelectedBond] = useState(null);
  const [bondId, setBondId] = useState(0);
  const [selectedAtoms, setSelectedAtoms] = useState([]);
  const [bondSmarts, setBondSmarts] = useState("-");
  const [nodeId, setNodeId] = useState(0);
  const [drawnSmarts, setDrawnSmarts] = useState("");
  const [drawSmiles, setDrawSmiles] = useState("CC(C)C(=O)O");

  useEffect(() => {
    const handleKeyDown = (event) => {
      const { key } = event;
      switch (key) {
        case "∫":
          setDrawMode("bond");
          break;
        case "å":
          setDrawMode("atom");
          break;
        case "Dead":
          setAtomSmarts("[N;H2;+0:" + currAtomMap + "]");
          break;
        case "ç":
          setAtomSmarts("[C;H2;+0:" + currAtomMap + "]");
          break;
        case "ø":
          setAtomSmarts("[O;H1;+0:" + currAtomMap + "]");
          break;
        case "ß":
          setAtomSmarts("[S;H1;+0:" + currAtomMap + "]");
          break;
        case "π":
          setAtomSmarts("[P;H1;+0:" + currAtomMap + "]");
          break;
        default:
          break;
      }
    };

    // if (tab === "drawer") {
    window.addEventListener("keydown", handleKeyDown);
    // } else {
    //   window.removeEventListener("keydown", handleKeyDown);
    // }
  }, []);

  function getSMARTSfromDrawing() {
    fetch("/api/get_smarts_from_drawing", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        drawerNodes: drawerNodes,
        drawerEdges: drawerEdges,
      }),
    })
      .then((response) => response.json())
      .then((data) => {
        setDrawnSmarts(data.smarts);
      });
  }

  useEffect(() => {
    getSMARTSfromDrawing();
  }, [drawerNodes, drawerEdges]);

  function renderDrawModeOptions() {
    if (drawMode === "atom") {
      return (
        <div className="sidebar-item">
          <div className="sidebar-item-label">atom smarts</div>
          <input
            className="sidebar-item-input"
            value={atomSmarts}
            onChange={(e) => setAtomSmarts(e.target.value)}
          ></input>
        </div>
      );
    }
    if (drawMode === "bond") {
      return (
        <div className="sidebar-item">
          <div className="sidebar-item-label">bond smarts</div>
          <input
            className="sidebar-item-input"
            value={bondSmarts}
            onChange={(e) => setBondSmarts(e.target.value)}
          ></input>
        </div>
      );
    }
  }

  function updateSelectedAtomSmarts(nodeId, atomSmarts) {
    let _drawerNodes = [...drawerNodes];
    for (let node of _drawerNodes) {
      if (node.id === nodeId) {
        node.atomSmarts = atomSmarts;
        break;
      }
    }

    setSelectedAtom({ ...selectedAtom, atomSmarts: atomSmarts });
    setDrawerNodes(_drawerNodes);
  }

  function updateSelectedBondSmarts(nodeId, bondSmarts) {
    let _drawerEdges = [...drawerEdges];
    for (let node of _drawerEdges) {
      if (node.id === nodeId) {
        node.bondSmarts = bondSmarts;
        break;
      }
    }

    setSelectedBond({ ...selectedBond, bondSmarts: bondSmarts });
    setDrawerEdges(_drawerEdges);
  }

  function deleteAtom(nodeId) {
    let _drawerNodes = drawerNodes.filter((node) => node.id !== nodeId);

    // delete bonds attached to atom
    let _drawerEdges = drawerEdges.filter(
      (edge) => edge.source !== nodeId && edge.target !== nodeId
    );

    setDrawerEdges(_drawerEdges);
    setDrawerNodes(_drawerNodes);
    setSelectedAtom(null);
  }

  function deleteBond(bondId) {
    let _drawerEdges = drawerEdges.filter((edge) => edge.id !== bondId);
    setDrawerEdges(_drawerEdges);
    setSelectedBond(null);
  }

  function clearDrawing() {
    setDrawerNodes([]);
    setDrawerEdges([]);
    setCurrAtomMap(0);
    setAtomSmarts("[C;H2;+0:0]");
  }

  function renderFromDrawSmiles() {
    fetch("/api/get_drawing_from_smiles", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        smiles: drawSmiles,
      }),
    })
      .then((response) => response.json())
      .then((data) => {
        if (data.error) {
          alert("invalid input");
          return;
        }

        let height = canvasDrawerRef.current.height;
        let width = canvasDrawerRef.current.width;

        let adjusted_nodes = data.output.nodes.map((node) => {
          let x = ((node.x + 4) * 0.5 * width * 0.5) / 3 + 500;
          let y = ((node.y + 4) * 0.5 * height * 0.5) / 3;

          console.log(x, y);
          return {
            ...node,
            x: x,
            y: y,
          };
        });

        setDrawerNodes(adjusted_nodes);
        setDrawerEdges(data.output.edges);
        setNodeId(data.output.nodes.length);
        setBondId(data.output.edges.length);
      });
  }

  function renderDrawModeSelect() {
    if (drawMode === "atom") {
      if (selectedAtom) {
        return (
          <div className="sidebar-item">
            <div className="sidebar-item-label">selected atom</div>
            <input
              className="sidebar-item-input"
              value={selectedAtom?.atomSmarts}
              onChange={(e) =>
                updateSelectedAtomSmarts(selectedAtom.id, e.target.value)
              }
            ></input>
            <button onClick={() => deleteAtom(selectedAtom.id)}>delete</button>
          </div>
        );
      }
    }
    if (drawMode === "bond") {
      if (selectedBond) {
        return (
          <div className="sidebar-item">
            <div className="sidebar-item-label">selected bond</div>
            <input
              className="sidebar-item-input"
              value={selectedBond?.bondSmarts}
              onChange={(e) =>
                updateSelectedBondSmarts(selectedBond.id, e.target.value)
              }
            ></input>
            <button onClick={() => deleteBond(selectedBond.id)}>delete</button>
          </div>
        );
      }
    }

    return <div className="sidebar-item"></div>;
  }

  function returnDrawingOptions() {
    return (
      <div className="sidebar2">
        <div className="sidebar-item">
          <div className="sidebar-item-label">draw mode</div>
          <select
            style={{ width: "60%" }}
            value={drawMode}
            onChange={(e) => setDrawMode(e.target.value)}
          >
            <option value="atom">atom</option>
            <option value="bond">bond</option>
          </select>
        </div>
        {renderDrawModeOptions()}
        {renderDrawModeSelect()}
        <div className="sidebar-item">{drawnSmarts}</div>
        <div className="sidebar-item">
          <input
            className="sidebar-item-input"
            value={drawSmiles}
            onChange={(e) => setDrawSmiles(e.target.value)}
          ></input>

          <button onClick={() => renderFromDrawSmiles()}>draw</button>
          <button onClick={() => clearDrawing()}>clear drawing</button>
        </div>
      </div>
    );
  }

  function returnCanvas() {
    return (
      <Drawer
        canvasRef={canvasDrawerRef}
        setDrawerNodes={setDrawerNodes}
        drawerNodes={drawerNodes}
        setDrawerEdges={setDrawerEdges}
        drawerEdges={drawerEdges}
        drawMode={drawMode}
        atomSmarts={atomSmarts}
        setAtomSmarts={setAtomSmarts}
        setSelectedAtom={setSelectedAtom}
        selectedAtom={selectedAtom}
        setSelectedBond={setSelectedBond}
        selectedBond={selectedBond}
        setNodeId={setNodeId}
        nodeId={nodeId}
        selectedAtoms={selectedAtoms}
        setSelectedAtoms={setSelectedAtoms}
        bondSmarts={bondSmarts}
        bondId={bondId}
        setBondId={setBondId}
        currAtomMap={currAtomMap}
        setCurrAtomMap={setCurrAtomMap}
      />
    );
  }

  return (
    <div className="drawerView">
      <div className="top">{returnCanvas()}</div>
      <div className="bottom">{returnDrawingOptions()}</div>
    </div>
  );
};

export default DrawerView;
