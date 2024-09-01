import React, { useState, useEffect, useRef } from "react";
import ReactDOM from "react-dom";
import "../static/main.css";
import DrawerView from "./DrawerView";
import CorpusView from "./CorpusView";
import SearchView from "./SearchView";

// ========================================
function Main() {
  const [tab, setTab] = useState("corpus");
  const [mechanisms, setMechanisms] = useState({});

  function returnView() {
    if (tab === "drawer") {
      return <DrawerView />;
    }
    if (tab === "corpus") {
      return <CorpusView mechanisms={mechanisms} />;
    }
    if (tab === "search") {
      return <SearchView mechanisms={mechanisms} />;
    }

  }

  function getAllTemplateSets() {
    fetch("/api/get_all_template_sets", {
      method: "GET",
    })
      .then((response) => response.json())
      .then((data) => {
        setMechanisms(data.mechanisms);
      });
  }

  useEffect(() => {
    getAllTemplateSets();
  }, []);

  return (
    <div className="container">
      <div className="header">
        <div className="navLeft">
          <div className="headerTitle">SMARTS Helper</div>
        </div>
        <div className="headerButtons">
        <div
            className={"nav-item" + (tab === "search" ? " active" : "")}
            onClick={() => setTab("search")}
          >
            search
          </div>
          <div
            className={"nav-item" + (tab === "drawer" ? " active" : "")}
            onClick={() => setTab("drawer")}
          >
            drawer
          </div>
          <div
            className={"nav-item" + (tab === "corpus" ? " active" : "")}
            onClick={() => setTab("corpus")}
          >
            corpus
          </div>
        </div>
      </div>
      <div className="view">{returnView()}</div>
    </div>
  );
}
ReactDOM.render(<Main />, document.getElementById("reactEntry"));
