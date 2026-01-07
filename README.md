# TableTools  
*A zero‑dependency Python library for explicit, predictable, fully user‑controlled tabular data workflows.*

TableTools is a pure‑Python toolkit for reading, transforming, analyzing, and visualizing structured data — all without external dependencies or hidden abstractions. It provides a complete, end‑to‑end environment for working with tabular data, built around two core components:

- **The TableTools API** — manages file I/O, directory handling, workflow utilities, and access to the toolbox ecosystem  
- **The Table object** — a dedicated structure for storing, cleaning, reshaping, analyzing, and visualizing tabular data  

Together, these components support a wide range of real‑world workflows, from data ingestion to final visualization.

---

## Features

### **Zero Dependencies**
TableTools runs anywhere Python runs — no installation, no environment conflicts, no external libraries.

### **Unified Top‑Level API**
The `TableTools` class provides:

- File input/output for many formats  
- Directory and file management  
- Workflow utilities (timers, progress counters, manual viewer)  
- Access to specialized toolboxes  

### **Powerful Table Object**
The `Table` object supports:

- Row/column access  
- Data cleaning and restructuring  
- Sampling and statistics  
- Sorting and filtering  
- Joins, splitting, comparison, pivoting, transposing  
- Expression‑based column calculations  
- Mapping/applying custom functions  
- Interpolation and filling  
- Built‑in HTML/SVG plotting  
- An interactive browser‑based spreadsheet viewer  

### **Extensive Format Support**
TableTools can read and write:

- CSV, TSV, and other delimited text  
- JSON  
- XML  
- YAML  
- HTML tables  
- SQL query results  
- Excel files  
- Fixed‑width text  
- dBase/xBase  
- R dput  
- Markdown  
- LaTeX  
- And other structured text formats  

All readers return a `Table` object.  
All writers accept one.

### **Modular Toolbox Ecosystem**
Specialized toolboxes extend TableTools into:

- List manipulation  
- Mathematical and statistical operations  
- Text and date processing  
- Unit conversion  
- Matrix algebra  
- Point‑based spatial operations  
- Synthetic data generation  
- Standalone plotting  

---

## Download & Setup

### **Download**
- Visit the TableTools GitHub page  
- Click **Code → Download ZIP**  
- Save the ZIP file to your preferred directory  

### **Setup**
- Unzip the downloaded file  
- Copy the `TableTools/` folder into the same directory as your Python script  
- Import and use:

```python
from TableTools.Table_Tools import TableTools

tt = TableTools()
```

No installation required.

---

## Key Objects

### **TableTools Object**
The top‑level API for:

- Reading/writing data  
- Managing input/output directories  
- Listing files and file sizes  
- Checking for file/directory existence  
- Creating/deleting files and directories  
- Controlling overwrite behavior  
- Running timers and progress counters  
- Opening the manual and help system  
- Accessing the toolbox ecosystem  

**Example:**

```python
from TableTools.Table_Tools import TableTools
tt = TableTools()
```

---

### **Table Object**
The primary data container in the library.

It stores tabular data in a structured row‑and‑column format and provides a large suite of methods for:

- Cleaning  
- Transforming  
- Reshaping  
- Analyzing  
- Visualizing  

#### **Accessing Data**

```python
row    = table[0]
column = table["column_name"]
value  = table[row_index, column_index]
```

This indexing is **read‑only**.  
All modifications are performed through Table methods.

#### **Capabilities**

- Metadata management  
- Row/column extraction  
- Sampling and statistics  
- Reshaping (pivot, transpose, split, join, compare)  
- Adding/removing/replacing rows and columns  
- Sorting  
- Filtering  
- Expression‑based column calculations  
- Mapping/applying custom functions  
- Interpolation and filling  
- HTML/SVG plotting  
- Interactive spreadsheet‑style viewer  

---

## Documentation

A full manual is included in:

```
TableTools/Docs/TableTools_manual.html
```

You can also open it programmatically:

```python
tt.view_manual()
```

---

## License

TableTools is released under the MIT License.  
See `LICENSE` for details.