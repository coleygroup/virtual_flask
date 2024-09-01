import React, { useEffect, useState, useRef } from "react";
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
  const [selectedIndices, setSelectedIndices] = useState([]);
  const [maxFreq, setMaxFreq] = useState(0);
  const [minFreq, setMinFreq] = useState(0);
  const [filtersApplied, setFiltersApplied] = useState([]);
  const [filterRemoved, setFilterRemoved] = useState([]);
  const [loading, setLoading] = useState(false);
  const [selectedDisplay, setSelectedDisplay] = useState(
    "min_delta_lambda_max"
  );

  function getHits() {
    setLoading(true);
    fetch("/api/get_hist", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({ selectedDisplay: selectedDisplay }),
    })
      .then((response) => response.json())
      .then((data) => {
        setHits(data.hist);
        setMaxFreq(data.max_freq);
        setLoading(false);
      });
  }

  function updateHits(hits) {
    fetch("/api/update_hist", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        ids: hits,
      }),
    })
      .then((response) => response.json())
      .then((data) => {
        setHits(data.hist);
        setMaxFreq(data.max_freq);
      });
  }

  useEffect(() => {
    getHits();
  }, []);

  const scrollToDiv = (divId) => {
    // Assuming divId is the id attribute of the div you want to scroll to
    const divElement = document.getElementById(divId);
    if (divElement) {
      divElement.scrollIntoView({ behavior: "smooth", block: "start" });
    }
  };

  function updateDrawerNodesAndEdges(
    path,
    mech_path,
    i,
    path_data,
    unmapped_product_smiles,
    id
  ) {
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
        unmapped_product_smiles: unmapped_product_smiles,
      }),
    })
      .then((response) => response.json())
      .then((data) => {
        setRenderNodes(data.nodes);
        setRenderEdges(data.edges);
        setRenderArrows(data.arrows);
        // console.log(data.intermediates)
        setIntermediates(data.intermediates);
        scrollToDiv(id);
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
        // console.log(data);
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

  function renderMechanism(row) {
    let output = [];
    console.log(row);
    console.log(pathData);
    for (let i = 0; i < pathData.length; i++) {
      if (intermediates[i + 1] == undefined) {
        // console.log("its happening")
        continue;
      }
      output.push(
        <div key={i} className="mechHolderInfo">
          <div className="mechHolder">
            <Renderer
              key={i}
              canvasRef={React.createRef()}
              drawerNodes={renderNodes[i]}
              drawerEdges={renderEdges[i]}
              drawerArrows={renderArrows[i]}
              idx={selectedProduct}
            />
            <img
              src={"data:image/png;base64," + intermediates[i + 1]}
              alt="missing smiles"
              style={{ height: "100%" }}
            />
          </div>
          <div className="mechInfo">
            {pathData[i].name}, {pathData[i].scope}, {pathData[i].description},{" "}
            {row[headers.indexOf("energies")][i].toString()}
            <br />
            <div className="smilesHolder">
              {row[headers.indexOf("mechanistic_path")][i + 1]}
            </div>
          </div>
        </div>
      );
    }
    return (
      <div className="mechBox">
        {/* <button onClick={() => downloadOutput(row)}>
          Download Assay Design
        </button> */}
        <div className="smilesHolder">
          {row[headers.indexOf("mapped_input_smiles")]}
        </div>
        {output}
      </div>
    );
  }

  function downloadOutput(sel) {
    let tit = "assay_design.xlsx";
    let worksheet = null;

    let headers2 = [
      "chemical_description",
      "concentration",
      "density",
      "smiles",
      "state_of_matter",
      "overhead",
      "factor",
      "reaction_concentration",
    ];
    // console.log(sel.mapped_input_smiles);

    let output = [];
    let idx1 = headers.indexOf("mapped_input_smiles");
    let sliced = sel[idx1].slice(1, -1);
    let split = sliced.split(",");

    for (let i = 0; i < split.length; i++) {
      let obj = {
        chemical_description: "component_" + i,
        concentration: 0.1,
        density: 1,
        smiles: split[i],
        state_of_matter: "liquid",
        overhead: 2,
        factor: "factor_" + i,
        reaction_concentration: 0.1,
      };
      output.push(obj);
    }

    // console.log(output);

    worksheet = XLSX.utils.json_to_sheet(output, {
      header: headers2,
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
    // console.log(output_smiles.length);

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

  const [results, setResults] = useState([]);

  const [filterString, setFilterString] = useState(
    'SELECT id FROM test_table5 WHERE NOT path_data @> \'[{"name": "Prins Reaction"}]\' AND NOT path_data @> \'[{"name": "Ene Reaction"}]\';'
  );

  const handleQuerySubmit = (event) => {
    event.preventDefault();
    setBufferIndex(0);
    setResults([]);
    setRenderNodes([]);
    setRenderEdges([]);
    setPathData([]);
    setRenderArrows([]);
    setIntermediates([]);
    setSelectedIndices([]);
    setHitBuffer([]);
    setSelectedBin(null);

    setLoading(true);
    fetch("/api/run_query", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        query: filterString,
        selectedDisplay: selectedDisplay,
      }),
    })
      .then((response) => response.json())
      .then((data) => {
        setHits(data.hist);
        setMaxFreq(data.max_freq);
        setFiltersApplied([...filtersApplied, filterString]);
        setFilterRemoved([...filterRemoved, data.removed]);
        setLoading(false);
        // setResults(data.found);
        // updateHits(data.found);
      });
  };

  const [headers, setHeaders] = useState([]);

  const renderAdditionalInfo = (rowIndex, row) => {
    return (
      <tr>
        <td key={rowIndex} colSpan={headers.length}>
          {renderMechanism(row)}
        </td>
      </tr>
    );
  };

  const TableComponent = ({ hitBuffer }) => {
    if (hitBuffer.length === 0) return null;
    if (headers.length === 0) return null;

    const skippers = [
      "mols",
      "morgan_fp",
      "mech_path",
      "path",
      "propagations",
      "mapped_product_smiles",
      "target_molecule",
      "unmapped_product_smiles",
      "path_data",
      "mapped_input_smiles",
      "unmapped_input_smiles",
      "overall_reaction",
      "mechanistic_path",
      "energies",
      "peakwavs_max_ensemble_uncal_var",
      "average_delta_g"
    ];
    return (
      <table style={{ width: "100%" }}>
        <thead>
          <tr>
            {headers.map((header, index) => {
              if (skippers.includes(header)) {
                return null;
              } else {
                return <th key={index}>{header}</th>;
              }
            })}
          </tr>
        </thead>
        <tbody>
          {hitBuffer.map((row, rowIndex) => (
            <React.Fragment key={rowIndex}>
              <tr
                id={`row-${rowIndex}`}
                key={rowIndex}
                onClick={() =>
                  updateDrawerNodesAndEdges(
                    row[headers.indexOf("mech_path")],
                    row[headers.indexOf("mechanistic_path")],
                    rowIndex,
                    row[headers.indexOf("path_data")],
                    row[headers.indexOf("unmapped_product_smiles")],
                    `row-${rowIndex}`
                  )
                }
              >
                {row.map((cell, colIndex) => {
                  const header = headers[colIndex];
                  if (header === "found_in_askcos_forward") {
                    if (cell === true) {
                      return <td key={colIndex}>true</td>;
                    } else if (cell === false) {
                      return <td key={colIndex}>false</td>;
                    }
                  }
                  if (skippers.includes(header)) {
                    return null;
                  }
                  if (header === "reaction_img") {
                    return (
                      <td key={colIndex}>
                        <img
                          src={`data:image/png;base64,${cell}`}
                          alt="Chemical reaction"
                          style={{ height: "80%" }}
                        />
                      </td>
                    );
                  }
                  return <td key={colIndex}>{cell}</td>;
                })}
              </tr>
              {selectedProduct === rowIndex &&
                renderAdditionalInfo(rowIndex, row)}
            </React.Fragment>
          ))}
        </tbody>
      </table>
    );
  };

  const [hitBuffer, setHitBuffer] = useState([]);
  const [bufferSize, setBufferSize] = useState(10);
  const [bufferIndex, setBufferIndex] = useState(0);
  const [selectedBin, setSelectedBin] = useState(null);

  function updateBuffer(indexes) {
    setLoading(true);
    fetch("/api/get_buffer", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        buffer: indexes,
      }),
    })
      .then((response) => response.json())
      .then((data) => {
        setHitBuffer(data.hits);
        setHeaders(data.headers);
        setLoading(false);
      });
  }

  function gridClickHandler(e, bin) {
    setBufferIndex(0);
    setSelectedBin(bin);
    setSelectedProduct(null);
    setRenderNodes([]);
    setRenderEdges([]);
    setPathData([]);
    setRenderArrows([]);
    setIntermediates([]);

    // console.log(bin.ids);

    if (bin.ids == null || bin.ids.length === 0) {
      setSelectedIndices([]);
    } else {
      setSelectedIndices(bin.ids);
      updateBuffer(bin.ids.slice(0, bufferSize));
    }

    // fetch("/api/get_bucket_data", {
    //   method: "POST",
    //   headers: {
    //     "Content-Type": "application/json",
    //   },
    //   body: JSON.stringify({
    //     bucket1: bin.bucket1,
    //     bucket2: bin.bucket2,
    //     ids: bin.ids,
    //   }),
    // })
    //   .then((response) => response.json())
    //   .then((data) => {
    //     if (data.hit_ids == null) {
    //       setSelectedIndices([]);
    //     } else {
    //       setSelectedIndices(data.hit_ids);
    //       updateBuffer(data.hit_ids.slice(0, bufferSize));
    //     }
    //   });
  }

  const getColor = (frequency) => {
    if (frequency === 0) return "rgba(200, 200, 200, 1)"; // Transparent white
    const alpha = Math.min(1, frequency / maxFreq + 0.50);
    return `rgba(31, 53, 105, ${alpha})`; // Red with variable opacity
  };

  const Grid = ({ data }) => {
    let num_buckets = Math.sqrt(data.length);
    console.log(num_buckets);
    const getBucketLabel = (index, axis) => {
      if (axis === "x") {
        return `${data[index].start1}-${data[index].end1}`;
      }
      return `Y: ${data[index].start2}-${data[index].end2}`;
    };
    return (
      <div
        style={{
          display: "flex",
          alignItems: "center",
          position: "relative",
        }}
      >
        <div
          style={{
            display: "flex",
            transform: "rotate(-90deg)",
            position: "absolute",
            left: "-60px",
          }}
        >
          max delta g
        </div>
        <div style={{ display: "flex", flexDirection: "column" }}>
          <div
            className="grid"
            style={{
              gridTemplateColumns: "repeat(" + num_buckets + ", 1fr)",
              width: "350px",
            }}
          >
            {data.map((bin, index) => (
              <div
                key={index}
                onClick={(e) => gridClickHandler(e, bin)}
                className="grid-cell"
                style={{
                  backgroundColor: getColor(bin.frequency),
                  gridColumnStart: bin.bucket1 + 1,
                  gridColumnEnd: bin.bucket1 + 2,
                  // gridRowStart: bin.bucket2,
                  // gridRowEnd: bin.bucket2 + 1,
                  gridRowStart: num_buckets - bin.bucket2,
                  gridRowEnd: num_buckets - (bin.bucket2 - 1),
                  // border: "1px solid #ccc",
                }}
                title={`Bucket: (${bin.start1}-${bin.end1}), (${bin.start2}-${bin.end2})\nFrequency: ${bin.frequency}`}
              />
            ))}
          </div>
          <div
            style={{
              display: "flex",
              width: "350px",
              justifyContent: "center",
            }}
          >
            {selectedDisplay}
          </div>
        </div>
      </div>
    );
  };

  function nextBuffer() {
    let working_array = selectedIndices;

    if (bufferIndex + bufferSize < working_array.length) {
      setSelectedProduct(null);
      setRenderNodes([]);
      setRenderEdges([]);
      setPathData([]);
      setRenderArrows([]);
      setIntermediates([]);

      setBufferIndex(bufferIndex + bufferSize);

      updateBuffer(
        working_array.slice(
          bufferIndex + bufferSize,
          bufferIndex + 2 * bufferSize
        )
      );
    }
  }

  function prevBuffer() {
    let working_array = selectedIndices;
    if (bufferIndex - bufferSize >= 0) {
      setSelectedProduct(null);
      setRenderNodes([]);
      setRenderEdges([]);
      setPathData([]);
      setRenderArrows([]);
      setIntermediates([]);

      setBufferIndex(bufferIndex - bufferSize);

      updateBuffer(working_array.slice(bufferIndex - bufferSize, bufferIndex));
    }
  }

  function returnTableButtons() {
    let bin1_info = "";
    let bin2_info = "";
    if (selectedBin != null) {
      bin1_info = `${selectedBin.start1} - ${selectedBin.end1}`;
      bin2_info = `${selectedBin.start2} - ${selectedBin.end2}`;
    }

    return (
      <div className="tableData">
        <div className="binInfoBox">
          <label>{bin1_info}</label>
          <label>{bin2_info}</label>
        </div>
        <div className="buttonBox2">
          <button onClick={nextBuffer}>next</button>
          <button onClick={prevBuffer}>prev</button>
        </div>
        {bufferIndex} - {bufferIndex + bufferSize} ({selectedIndices.length}{" "}
        hits)
      </div>
    );
  }

  function clearQuery() {
    setResults([]);
    setBufferIndex(0);
    updateBuffer(selectedIndices.slice(0, bufferSize));
  }

  function renderFilters() {
    let filts = [];
    console.log(filtersApplied);
    for (let i = 0; i < filtersApplied.length; i++) {
      filts.push(
        <div key={i}>
          {filtersApplied[i]}, {filterRemoved[i]}
        </div>
      );
    }
    return <div className="filtersApplied">{filts}</div>;
  }

  function renderLoadingBox() {
    if (loading) {
      return (
        <div className="loadingBox">
          <div className="loadingText">Loading...</div>
        </div>
      );
    }
    return;
  }


  useEffect(() => {
    getHits();
  }, [selectedDisplay]);

  function handleSelectChange(v){
    setSelectedDisplay(v);
  }

  function renderDisplayOptions() {
    return (
      <select
        value={selectedDisplay}
        onChange={(e) => handleSelectChange(e.target.value)}
      >
        <option value="min_delta_lambda_max">min_delta_lambda_max</option>
        <option value="network_size">network_size</option>
      </select>
    );
  }

  return (
    <div className="container">
      {renderLoadingBox()}
      {renderDisplayOptions()}
      <Grid data={hits} />

      <div className="queryBox">
        <form onSubmit={handleQuerySubmit} className="formClass">
          <div className="filterBox">
            <textarea
              type="text"
              value={filterString}
              onChange={(e) => setFilterString(e.target.value)}
            />
          </div>
          <div className="buttonBox">
            <button type="submit">Run Query</button>
          </div>
        </form>
        {renderFilters()}
      </div>

      <div className="resultsBox">
        {returnTableButtons()}
        <TableComponent hitBuffer={hitBuffer} />
      </div>
    </div>
  );
}

// ========================================

ReactDOM.render(<Main />, document.getElementById("reactEntry"));
