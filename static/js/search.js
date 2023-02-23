// Grabbing the inputs
let Meta = document.getElementById("search_meta_select")
let form = document.getElementById("search_form");
let smileSearch = document.getElementById("search_SMILES");
let smileSearchHBD = document.getElementById("search_SMILES_HBD");
let smileSearchHBA = document.getElementById("search_SMILES_HBA")
let mwSearch = document.getElementById("search_MW");
let hbdMwSearch = document.getElementById("search_HBD_MW")
let hbaMwSearch = document.getElementById("search_HBA_MW")




var selectedMeta = Meta.value
// Choosing which type of query you want to run, will show or hide that particular form
const setMeta = () => {
  selectedMeta = Meta.value
  let moleculeInputs = document.getElementsByClassName("search-molecule")
  let desInputs = document.getElementsByClassName("search-des")
  if (selectedMeta === "search-query-des") {
    moleculeInputs[0].style.display = "none"
    moleculeInputs[0].disabled = true;
    desInputs[0].style.display = "block"
    desInputs[0].disabled = false;


  }
  else {
    moleculeInputs[0].style.display = "block"
    moleculeInputs[0].disabled = false;
    desInputs[0].style.display = "none"
    desInputs[0].disabled = true;

  }
}
// using the selected radio button, wil display/hide the appropriate SMILES input bar
const displaySmileSearch = (event) => {
  if (event.target.value == "Yes") {
    if (event.target.name === "smiles_search") {
      smileSearch.style.display = "flex"
    }
    else if (event.target.name === "smiles_search_hbd") {
      smileSearchHBD.style.display = "flex"
    }
    else {
      smileSearchHBA.style.display = "flex"
    }

  }
  else {
    if (event.target.name === "smiles_search") {
      smileSearch.style.display = "none"
    }
    else if (event.target.name === "smiles_search_hbd") {
      smileSearchHBD.style.display = "none"
    }
    else {
      smileSearchHBA.style.display = "none"
    }
  }
}

// using the selected radio button, wil display/hide the appropriate MW input bar
const displayMwSearch = (event) => {
  console.log(event.target.name)
  if (event.target.value == "Yes") {

    if (event.target.name === "mw_search") {
      mwSearch.style.display = "flex"
    }
    else if (event.target.name === "hbd_mw_search") {
      hbdMwSearch.style.display = "flex"
    }
    else {
      hbaMwSearch.style.display = "flex"

    }
  }
  else {
    if (event.target.name === "mw_search") {
      mwSearch.style.display = "none"
    }
    else if (event.target.name === "hbd_mw_search") {
      hbdMwSearch.style.display = "none"
    }
    else {
      hbaMwSearch.style.display = "none"
    }
  }

}

//generates table head for search output
const generateTableHead = (table, data) => {
  let thead = table.createTHead();
  let row = thead.insertRow();
  for (let key of data) {
    let th = document.createElement("th");
    let text = document.createTextNode(key);
    th.appendChild(text);
    row.appendChild(th);
  }
}

// generates table for search output
const generateTable = (div, sdb) => {
  let table = div.appendChild(document.createElement("table"))
  table.setAttribute("id", "search-table-table")

  if (sdb["data"] === undefined || sdb["data"].length == 0) {
    let p = document.createElement('p')
    p.innerHTML = "No Such Data Found"
    div.appendChild(p)
    return
  }
  generateTableHead(table, sdb["columns"])
  let tbody = table.createTBody()
  for (let element of sdb["data"]) {
    let row = tbody.insertRow();
    for (key in element) {
      let cell = row.insertCell();
      let text = document.createTextNode(element[key]);
      cell.appendChild(text);
    }
  }

  $('#search-table-table').DataTable(
    { "lengthMenu": [5, 10, 15] }
  );

}


// Before form submission remove keys that have empty stings
document.getElementById("search_form").addEventListener("submit", function (event) {
  event.preventDefault();

  let toRemove = [];
  for (let i = 0; i < this.elements.length; i++) {

    if (this.elements[i].value === "" || this.elements[i].value === "No") {
      toRemove.push(this.elements[i])
    }


  }
  for (i of toRemove) {
    i.remove()
  }
  console.log("hey")

  // Submit the form
  this.submit();
});

// Initially run setMeta to hide the Des form
setMeta();

// Initially the MW inputs and smiles are hidden
smileSearch.style.display = "none"
smileSearchHBA.style.display = "none"
smileSearchHBD.style.display = "none"
mwSearch.style.display = "none"
hbdMwSearch.style.display = "none"
hbaMwSearch.style.display = "none"