import React, { useState, useEffect, useRef } from "react";
import Renderer from "./renderer";

const CorpusView = (mechanisms) => {
  const [selectedMechanism, setSelectedMechanism] = useState(null);
  const [selectedScope, setSelectedScope] = useState(null);
  const [renderNodes, setRenderNodes] = useState([]);
  const [renderEdges, setRenderEdges] = useState([]);
  const [renderArrows, setRenderArrows] = useState([]);

  const [testInputs, setTestInputs] = useState([]);
  const [testOutputs, setTestOutputs] = useState([]);

  function updateDrawerNodesAndEdges(mechanism, scope) {
    fetch("/api/get_drawer_nodes_and_edges", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        mechanism: mechanisms.mechanisms[mechanism].mechanism,
        scope: scope,
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

  function changeSelectedMechanism(mechanism) {
    // console.log(mechanism);
    setSelectedMechanism(mechanism);
    setSelectedScope(null);
    updateDrawerNodesAndEdges(mechanism, null);
    let _testInputs = [];
    for (
      let i = 0;
      i < mechanisms.mechanisms[mechanism].mechanism.length;
      i++
    ) {
      _testInputs.push("");
    }
    setTestInputs(_testInputs);
    let _testOutputs = [];
    for (
      let i = 0;
      i < mechanisms.mechanisms[mechanism].mechanism.length;
      i++
    ) {
      _testOutputs.push("");
    }
    setTestOutputs(_testOutputs);
  }

  function changeSelectedScope(scope) {
    setSelectedScope(scope);
    updateDrawerNodesAndEdges(selectedMechanism, scope);
  }

  function renderMechanisms() {
    let mechanismList = [];
    for (let mechanism in mechanisms.mechanisms) {
      let col = "black";
      if (selectedMechanism === mechanism) {
        col = "red";
      }
      mechanismList.push(
        <div className="mechanism" key={mechanism}>
          <div
            style={{ color: col }}
            key={mechanism}
            className="mechanismTitle"
            onClick={() => changeSelectedMechanism(mechanism)}
          >
            {mechanism}
          </div>
        </div>
      );
    }
    return mechanismList;
  }

  function queryPattern(step) {
    fetch("/api/query_pattern", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        mechanism: mechanisms.mechanisms[selectedMechanism].mechanism[step],
        scope: selectedScope,
      }),
    })
      .then((response) => response.json())
      .then((data) => {
        // set input to output of api
        setTestInputs((testInputs) => {
          let newTestInputs = [...testInputs];
          newTestInputs[step] = data.output;
          return newTestInputs;
        });
      });
  }

  function runTestInput(step, input, intra) {
    fetch("/api/run_test_input", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        mechanism: mechanisms.mechanisms[selectedMechanism].mechanism[step],
        scope: selectedScope,
        input: input,
        intra: intra,
      }),
    })
      .then((response) => response.json())
      .then((data) => {
        console.log(data.output, step);
        setTestOutputs((testOutputs) => {
          let newTestOutputs = [...testOutputs];
          newTestOutputs[step] = data.output;
          return newTestOutputs;
        });
      });
  }

  function renderSelectedMechanism() {
    if (selectedMechanism) {
      const scope = [];
      for (
        let i = 0;
        i < mechanisms.mechanisms[selectedMechanism].scope.length;
        i++
      ) {
        let content = [];
        for (let j in mechanisms.mechanisms[selectedMechanism].scope[i]) {
          content.push(
            <div key={j} className="mechanismScopeItemContent">
              {j}: {mechanisms.mechanisms[selectedMechanism].scope[i][j]}
            </div>
          );
        }

        let col = "black";
        if (
          selectedScope === mechanisms.mechanisms[selectedMechanism].scope[i]
        ) {
          col = "red";
        }

        scope.push(
          <div
            className="mechanismScopeContent"
            style={{ color: col }}
            key={i}
            onClick={() =>
              changeSelectedScope(
                mechanisms.mechanisms[selectedMechanism].scope[i]
              )
            }
          >
            {content}
          </div>
        );
      }

      const steps = [];
      for (
        let i = 0;
        i < mechanisms.mechanisms[selectedMechanism].mechanism.length;
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
                  mechanisms.mechanisms[selectedMechanism].mechanism[i][
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
                mechanisms.mechanisms[selectedMechanism].mechanism[i][
                  "description"
                ]
              }
              )
              <input
                type="text"
                value={testInputs[i]}
                onChange={(e) => {
                  setTestInputs((testInputs) => {
                    let newTestInputs = [...testInputs];
                    newTestInputs[i] = e.target.value;
                    return newTestInputs;
                  });
                }}
              />
              <button onClick={() => queryPattern(i)}>Query</button>
              {/* <button onClick={() => runTestInput(i, testInputs[i], false)}>
                Test
              </button> */}
              <button onClick={() => runTestInput(i, testInputs[i], true)}>
                Test
              </button>
              <br />
              {testOutputs[i]}
            </div>
          </div>
        );
      }

      return (
        <div className="selectedMechanism">
          <div className="mechanismName">{selectedMechanism}</div>
          <div className="mechanismScopes">{scope}</div>
          <div className="mechanismSteps">{steps}</div>
        </div>
      );
    }
  }

  return (
    <div className="corpusView">
      <div className="mechanismList">{renderMechanisms()}</div>
      <div className="selectedMechanismViewer">{renderSelectedMechanism()}</div>
    </div>
  );
};

export default CorpusView;
