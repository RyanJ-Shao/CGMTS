# CGM Time Series (CGMTS) User Guide
## 1.	Packages which CGMTS depends on.
```
dplyr
plotly
imputeTS
forecast
TSA
```

## 2.	Install “orca”
CGMTS package needs orca to export pdf file of CGM plot. Orca is an Electron app that generates images and reports of Plotly things like plotly.js graphs, dash apps, dashboards from the command line. The installation guide of orca can be seen in https://github.com/plotly/orca.
Installation guide of orca in MacOS:
•	Unzip the mac-release.zip file.
•	Double-click on the orca-X.Y.Z.dmg file. This will open an installation window.
•	Drag the orca icon into the Applications folder.
•	Open finder and navigate to the Applications/ folder.
•	Right-click on the orca icon and select Open from the context menu.
•	A password dialog will appear asking for permission to add orca to your system PATH.
•	Enter your password and click OK.
•	This should open an Installation Succeeded window.
•	Open a new terminal and verify that the orca executable is available on your PATH.
Installation guide of orca in Windows:
•	Extract the windows-release.zip file.
•	In the release folder, double-click on orca Setup X.Y.Z, this will create an orca icon on your Desktop.
•	Right-click on the orca icon and select Properties from the context menu.
•	From the Shortcut tab, copy the directory in the Start in field.
•	Add this Start in directory to your system PATH (see below).
•	Open a new Command Prompt and verify that the orca executable is available on your PATH.
