// function to get tables of selected database
const getTableNames = async () => {
    let selectElement = document.getElementsByName("table_selected")[0]
    const selectedVal = document.getElementsByName("db_selected")[0].value

    fetch("/fetchTableNames", {
        method: "POST",
        body: JSON.stringify({ "db_name": selectedVal }),
        headers: {
            "Content-Type": "application/json",
        },
    })
        .then((res) => res.json())
        .then((res) => {
            let eleList = [];
            res["cols"].forEach((item) => {
                let opt = document.createElement('option')
                opt.value = item;
                opt.innerHTML = item;
                eleList.push(opt)

            })
            selectElement.replaceChildren(...eleList)
        })
}

// calling the function to fetch values for the first table

getTableNames()


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
    div.replaceChildren()
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
        { "lengthMenu": [5, 10, 15], scrollX: true, }
    );

}
const fetchTable = (e) => {
    e.preventDefault();
    let selectedVal = document.getElementsByName("table_selected")[0].value
    let dbName = document.getElementsByName("db_selected")[0].value
    let dispTableElement = document.getElementById("view-table")
    fetch("/fetchTable", {
        method: "POST",
        body: JSON.stringify({ "table_name": selectedVal, "db_name": dbName }),
        headers: {
            "Content-Type": "application/json",
        },
    })
        .then((res) => res.json())
        .then((res) => {
            generateTable(dispTableElement, res)
        })
}