import React, { useEffect, useState } from "react";
import ReactDOM, { render } from "react-dom";
import "../static/main.css";
import Renderer from "./renderer"; // Adjust the import path as necessary
import ClickableCanvas from "./ClickableCanvas"; // Adjust the import path as necessary
import * as XLSX from "xlsx";

function Main() {
  const [hits, setHits] = useState([]);
  const [selected, setSelected] = useState(null);
  const [renderNodes, setRenderNodes] = useState([]);
  const [renderEdges, setRenderEdges] = useState([]);
  const [renderArrows, setRenderArrows] = useState([]);
  const [intermediates, setIntermediates] = useState([]);
  const [selectedProduct, setSelectedProduct] = useState(null);
  const [pathData, setPathData] = useState([]);
  const [inputs, setInputs] = useState([]);

  function getHits() {
    fetch("/api/get_hits", {
      method: "GET",
    })
      .then((response) => response.json())
      .then((data) => {
        setHits(data.hits);
        setInputs(data.inputs);
        let smiles = [];
        for (let i = 0; i < data.inputs.length; i++) {
          smiles.push(data.inputs[i].smiles);
        }

        setRemaningInputs(smiles);
      });
  }

  useEffect(() => {
    getHits();
  }, []);

  function updateDrawerNodesAndEdges(path, mech_path, i, path_data) {
    console.log(i, path_data);
    setSelectedProduct(i);
    setPathData(path_data);
    fetch("/api/get_drawer_nodes_and_edges", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        mechanism: path,
        intermediates: mech_path,
      }),
    })
      .then((response) => response.json())
      .then((data) => {
        setRenderNodes(data.nodes);
        setRenderEdges(data.edges);
        setRenderArrows(data.arrows);
        setIntermediates(data.intermediates);
      });
  }

  function getFromS3(smiles) {
    fetch("/api/get_from_s3", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        smiles: smiles,
      }),
    })
      .then((response) => response.json())
      .then((data) => {
        if (data == null) {
          setSelected([]);
        } else {
          // console.log(data.data);
          setSelected(data.data);
        }
      });
  }

  const [remainingInputs, setRemaningInputs] = useState([]);

  function updateList(smiles) {
    fetch("/api/update_list", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        smiles: smiles,
      }),
    })
      .then((response) => response.json())
      .then((data) => {
        console.log(data);
        setRemaningInputs(data);
      });
  }

  function updateSelected(e) {
    getFromS3(e);
    setSelectedProduct(null);
    setRenderNodes([]);
    setRenderEdges([]);
    setPathData([]);
    setRenderArrows([]);
    setIntermediates([]);
    setSelected(e);
  }

  function renderSelected() {
    if (!selected) {
      return <p>No selection</p>;
    }
    if (selected.length == 0) {
      return <p>No hits</p>;
    }
    let outputs = [];
    let sele = "";
    if ("data" in selected) {
      for (let i = 0; i < selected.data.length; ) {
        let item = selected.data[i];
        if (i === selectedProduct) {
          sele = " selected";
        } else {
          sele = "";
        }
        outputs.push(
          <div key={i} className={"productHolder" + sele}>
            <div className="productInfo">
              {item.unmapped_product_smiles}{" "}
              <button
                onClick={() =>
                  updateDrawerNodesAndEdges(
                    item.path,
                    item.mech_path,
                    i - 1,
                    item.path_data
                  )
                }
              >
                check mech.
              </button>
            </div>
            <img
              src={"data:image/png;base64," + item.img}
              alt="missing smiles"
              style={{ height: "80%" }}
            />
            {/* {item.path} */}
          </div>
        );
        i = i + 1;
      }
    }
    return (
      <div className="selectedHolder">
        <div className="outputsHolder">{outputs}</div>
      </div>
    );
  }

  function renderMechanism() {
    let output = [];
    if (selectedProduct != null) {
      intermediates.push(selected.data[selectedProduct].img);
    }
    console.log(intermediates.length, pathData.length);
    for (let i = 0; i < pathData.length; i++) {
      output.push(
        <div className="mechHolderInfo">
          <div className="mechHolder">
            <Renderer
              key={i}
              canvasRef={React.createRef()}
              drawerNodes={renderNodes[i]}
              drawerEdges={renderEdges[i]}
              drawerArrows={renderArrows[i]}
            />
            <img
              src={"data:image/png;base64," + intermediates[i + 1]}
              alt="missing smiles"
              // style={{ width: "30%" }}
            />
          </div>
          <div className="mechInfo">
            {pathData[i].name}, {pathData[i].scope}, {pathData[i].description}
          </div>
        </div>
      );
    }
    return <div className="mechBox">{output}</div>;
  }

  function downloadOutput() {
    let tit = "assay_design.xlsx";
    let worksheet = null;

    let headers = [
      "chemical_description",
      "concentration",
      "density",
      "smiles",
      "state_of_matter",
      "overhead",
      "factor",
      "reaction_concentration",
    ];
    console.log(selected.mapped_input_smiles);

    let output = [];
    for (let i = 0; i < selected.mapped_input_smiles.length; i++) {
      let obj = {
        chemical_description: "component_" + i,
        concentration: 0.1,
        density: 1,
        smiles: selected.mapped_input_smiles[i],
        state_of_matter: "liquid",
        overhead: 2,
        factor: "factor_" + i,
        reaction_concentration: 0.1,
      };
      output.push(obj);
    }

    console.log(output);

    worksheet = XLSX.utils.json_to_sheet(output, {
      header: headers,
    });
    const workbook = XLSX.utils.book_new();
    XLSX.utils.book_append_sheet(workbook, worksheet, "Sheet1");

    const workbookOut = XLSX.write(workbook, {
      bookType: "xlsx",
      type: "array",
    });

    const blob = new Blob([workbookOut], { type: "application/octet-stream" });

    const link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = tit;

    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  }

  function downloadMasses() {
    let headers = [
      "product_smiles",
      "input_smiles",
      "m/z",
      "M+1", // +H
      "M-1", // -H
      "M+18", // +NH4
      "M+23", // +Na
      "M+39", // +K
      "M+33", // +CH3OH+H
      "M+42", // +ACN+H
      "M+64", // +ACN+Na
      "M+79", // +DMSO+H
      "M+83", // +2ACN+H
    ];

    let output_smiles = [];
    if ("data" in selected) {
      for (let i = 0; i < selected.data.length; i++) {
        output_smiles.push(selected.data[i].mapped_product_smiles);
      }
    }
    console.log(output_smiles.length);

    fetch("/api/get_product_masses", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        smiles: output_smiles,
      }),
    })
      .then((response) => response.json())
      .then((data) => {
        let tit = "product_masses.xlsx";
        let worksheet = null;

        let input_smiles = [];
        for (let i = 0; i < selectedIndices.length; i++) {
          input_smiles.push(inputs[selectedIndices[i]].smiles);
        }

        let output = [];
        for (let i = 0; i < data.masses.length; i++) {
          let obj = {
            product_smiles: output_smiles[i],
            input_smiles: input_smiles.join("."),
            "m/z": data.masses[i],
            "M+1": data.masses[i] + 1,
            "M-1": data.masses[i] - 1,
            "M+18": data.masses[i] + 18,
            "M+23": data.masses[i] + 23,
            "M+39": data.masses[i] + 39,
            "M+33": data.masses[i] + 33,
            "M+42": data.masses[i] + 42,
            "M+64": data.masses[i] + 64,
            "M+79": data.masses[i] + 79,
            "M+83": data.masses[i] + 83,
          };
          output.push(obj);
        }
        worksheet = XLSX.utils.json_to_sheet(output, {
          header: headers,
        });
        const workbook = XLSX.utils.book_new();
        XLSX.utils.book_append_sheet(workbook, worksheet, "Sheet1");

        const workbookOut = XLSX.write(workbook, {
          bookType: "xlsx",
          type: "array",
        });

        const blob = new Blob([workbookOut], {
          type: "application/octet-stream",
        });

        const link = document.createElement("a");
        link.href = URL.createObjectURL(blob);
        link.download = tit;

        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
      });
  }

  // function renderInputGrid() {
  //   let grid = [];
  //   console.log(inputs)
  //   for (let i = 0; i < inputs.length; i++) {
  //     grid.push(
  //       <div>
  //         <img
  //           src={"data:image/png;base64," + inputs[i].img}
  //           alt="missing smiles"
  //           style={{ width: "30%" }}
  //         />
  //       </div>
  //     );
  //   }
  //   return <div className="inputGrid"> {grid}</div>;
  // }

  const [selectedIndices, setSelectedIndices] = useState([]);

  // Click handler to manage selected images
  const handleImageClick = (index, new_smiles) => {
    setSelectedProduct(null);
    setRenderNodes([]);
    setRenderEdges([]);
    setPathData([]);
    setRenderArrows([]);
    setIntermediates([]);

    if (selectedIndices.includes(index)) {
      // If already selected, remove from the selected array
      let smiles = [];
      for (let i = 0; i < selectedIndices.length; i++) {
        if (inputs[selectedIndices[i]].smiles == new_smiles) {
          continue;
        }
        smiles.push(inputs[selectedIndices[i]].smiles);
      }
      updateList(smiles);
      setSelectedIndices(selectedIndices.filter((i) => i !== index));
    } else {
      if (selectedIndices.length > 2) {
      } else if (selectedIndices.length == 2) {
        let smiles = [];
        for (let i = 0; i < selectedIndices.length; i++) {
          smiles.push(inputs[selectedIndices[i]].smiles);
        }
        smiles.push(new_smiles);

        getFromS3(smiles);
        setSelectedIndices([...selectedIndices, index]);
      } else {
        let smiles = [];
        for (let i = 0; i < selectedIndices.length; i++) {
          smiles.push(inputs[selectedIndices[i]].smiles);
        }
        smiles.push(new_smiles);
        updateList(smiles);
        setSelectedIndices([...selectedIndices, index]);
      }
      // Add to the selected array
    }
  };

  function renderInputGrid() {
    const grid = inputs.map((input, index) => {
      if (!remainingInputs.includes(input.smiles)) {
        return null;
      }

      return (
        <div
          key={index}
          onClick={() => handleImageClick(index, input.smiles)} // Add click handler
          style={{
            border: "1px solid red",
            // width: '19%',
            border: selectedIndices.includes(index)
              ? "2px solid red"
              : "2px solid black", // Conditional border style
            overflow: "hidden",
            cursor: "pointer",
          }}
        >
          <img
            src={`data:image/png;base64,${input.img}`}
            alt="Displayed Image"
            style={{
              width: "100%",
            }}
          />
        </div>
      );
    });

    return (
      <div
        className="inputGrid"
        style={{
          display: "flex",
          flexWrap: "wrap",
          width: "100%",
          height: "100%",
          overflowY: "auto",
        }}
      >
        {grid}
      </div>
    );
  }

  function renderSelected2() {
    const selectedImages = selectedIndices.map((index) => (
      <div>
        <img
          src={`data:image/png;base64,${inputs[index].img}`}
          alt="Selected Image"
          style={{
            width: "100%",
          }}
        />
      </div>
    ));

    return (
      <div
        className="selectedImages"
        style={{
          display: "flex",
          flexWrap: "wrap",
          width: "100%",
          height: "90%",
          overflowY: "auto",
        }}
      >
        {selectedImages}
      </div>
    );
  }

  function returnSmiles() {
    let smiles = [];
    for (let i = 0; i < selectedIndices.length; i++) {
      smiles.push(inputs[selectedIndices[i]].smiles);
      if (i != selectedIndices.length - 1) {
        smiles.push(".");
      }
    }

    return smiles;
  }

  
  // if ("data" in selected) {
  //   let sel_len = selected.data.length;
  // }

  return (
    <div className="container">
      <div className="left">
        <div className="selectedHolder">
          {renderSelected2()}
          {returnSmiles()}
          <button onClick={() => downloadOutput()}>download assay input</button>
          <button onClick={() => downloadMasses()}>
            download assay output masses
          </button>
        </div>
        <div className="canvasHolder">
          {/* <ClickableCanvas
            data={hits}
            setSelected={updateSelected}
            selected={selected}
          /> */}
          {renderInputGrid()}
        </div>
      </div>
      <div className="right">
        {renderSelected()}
        {renderMechanism()}
      </div>
    </div>
  );
}

// ========================================

ReactDOM.render(<Main />, document.getElementById("reactEntry"));
