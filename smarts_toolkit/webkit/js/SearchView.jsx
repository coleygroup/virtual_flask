import React, { useState, useEffect, useRef } from "react";
import Renderer from "./renderer";

const SearchView = (mechanisms) => {
  const [selectedMechanism, setSelectedMechanism] = useState(null);
  const [selectedScope, setSelectedScope] = useState(null);
  const [renderNodes, setRenderNodes] = useState([]);
  const [renderEdges, setRenderEdges] = useState([]);
  const [renderArrows, setRenderArrows] = useState([]);

  const [testInputs, setTestInputs] = useState([]);
  const [testOutputs, setTestOutputs] = useState([]);

  const [searchBar, setSearchBar] = useState("N[C@H](C(O)=O)CC1=CNC2=C1C=CC=C2.O=CC3=CC=CC=C3>>O=C(C(/[NH+]=C/C4=CC=CC=C4)CC5=CNC6=C5C=CC=C6)O");
  const [hitMechs, setHitMechs] = useState([])

  function updateDrawerNodesAndEdges(hits) {
    fetch("/api/get_drawer_nodes_and_edges", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        mechanism: hits,
        scope: null,
      }),
    })
      .then((response) => response.json())
      .then((data) => {
        setRenderNodes(data.nodes);
        setRenderEdges(data.edges);
        // console.log(data.arrows);
        setRenderArrows(data.arrows);
      });
  }


  function renderSelectedMechanism() {
    if (hitMechs) {
      const steps = [];
      for (
        let i = 0;
        i < hitMechs.length;
        i++
      ) {
        if (renderNodes[i] === undefined || renderEdges[i] === undefined) {
          continue;
        }
        steps.push(
          <div key={i} className="mechanismStepContent">
            <div className="rendererHolder">
              <div className="rendererMechanismText">
                {
                  hitMechs[i][
                    "template"
                  ]
                }
              </div>
              <Renderer
                key={i}
                canvasRef={React.createRef()}
                drawerNodes={renderNodes[i]}
                drawerEdges={renderEdges[i]}
                drawerArrows={renderArrows[i]}
                scope={selectedScope}
              />
            </div>
            <div className="testInput">
              (
              {
                hitMechs[i][
                  "description"
                ]
              }
              )
            </div>
          </div>
        );
      }

      if (steps.length === 0) {
        return <div className="selectedMechanism">No mechanism found</div>;
      }

      return (
        <div className="selectedMechanism">
          <div className="mechanismName">{selectedMechanism}</div>
          <div className="mechanismSteps">{steps}</div>
        </div>
      );
    }
  }


  function searchMechanisms(){
    fetch("/api/search_mechanisms", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        search: searchBar,
        mechanisms: mechanisms,
      }),
    })
      .then((response) => response.json())
      .then((data) => {
        console.log(data);
        setHitMechs(data.hits);
        updateDrawerNodesAndEdges(data.hits);
      });
  }

  function renderSearchBar(){
    return (
      <div className="searchBar">
        <input
          type="text"
          value={searchBar}
          style={{ width: "80%" }}
          onChange={(e) => {
            setSearchBar(e.target.value);
          }}
        />
        <button onClick={() => searchMechanisms()}>Search</button>
      </div>
    );
  }

  return (
    <div className="searchView">
      {renderSearchBar()}
      {/* <div className="mechanismList">{renderMechanisms()}</div> */}
      <div className="selectedMechanismViewer search">{renderSelectedMechanism()}</div>
    </div>
  );
};

export default SearchView;
