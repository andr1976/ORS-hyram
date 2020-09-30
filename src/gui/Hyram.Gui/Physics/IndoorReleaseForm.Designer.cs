﻿using System;
using System.ComponentModel;
using System.Drawing;
using System.Windows.Forms;

namespace SandiaNationalLaboratories.Hyram {
	partial class IndoorReleaseForm {

		#region Component Designer generated code
        private void InitializeComponent() {
            System.ComponentModel.ComponentResourceManager resources = new System.ComponentModel.ComponentResourceManager(typeof(IndoorReleaseForm));
            System.Windows.Forms.DataGridViewCellStyle dataGridViewCellStyle2 = new System.Windows.Forms.DataGridViewCellStyle();
            this.IOTabs = new System.Windows.Forms.TabControl();
            this.InputTab = new System.Windows.Forms.TabPage();
            this.PhaseSelection = new System.Windows.Forms.ComboBox();
            this.PhaseLabel = new System.Windows.Forms.Label();
            this.OverpressureSpinner = new System.Windows.Forms.PictureBox();
            this.InputWarning = new System.Windows.Forms.Label();
            this.InputTabs = new System.Windows.Forms.TabControl();
            this.IndoorReleaseTab = new System.Windows.Forms.TabPage();
            this.GeometryPicture = new System.Windows.Forms.PictureBox();
            this.InputGrid = new System.Windows.Forms.DataGridView();
            this.OutputOptionsTab = new System.Windows.Forms.TabPage();
            this.MaxTimeInput = new System.Windows.Forms.TextBox();
            this.MaxTimeLabel = new System.Windows.Forms.Label();
            this.PlottingOptionsGroupBox = new System.Windows.Forms.GroupBox();
            this.PressuresPerTimeCheckbox = new System.Windows.Forms.CheckBox();
            this.PressuresPerTimeGroupBox = new System.Windows.Forms.GroupBox();
            this.PressuresPerTimeGrid = new System.Windows.Forms.DataGridView();
            this.MarkTimeGridCol = new System.Windows.Forms.DataGridViewTextBoxColumn();
            this.MarkPressureGridCol = new System.Windows.Forms.DataGridViewTextBoxColumn();
            this.PressureLinesGroupBox = new System.Windows.Forms.GroupBox();
            this.PressureLinesGrid = new System.Windows.Forms.DataGridView();
            this.PressureLineCol = new System.Windows.Forms.DataGridViewTextBoxColumn();
            this.PressureLinesCheckbox = new System.Windows.Forms.CheckBox();
            this.PlotTimesInput = new System.Windows.Forms.TextBox();
            this.PlotTimesLabel = new System.Windows.Forms.Label();
            this.ExecuteBtn = new System.Windows.Forms.Button();
            this.outputTab = new System.Windows.Forms.TabPage();
            this.splitContainer1 = new System.Windows.Forms.SplitContainer();
            this.tcOutput = new System.Windows.Forms.TabControl();
            this.tpPressure = new System.Windows.Forms.TabPage();
            this.label2 = new System.Windows.Forms.Label();
            this.tpFlammableMass = new System.Windows.Forms.TabPage();
            this.tpLayer = new System.Windows.Forms.TabPage();
            this.tpTrajectory = new System.Windows.Forms.TabPage();
            this.tpData = new System.Windows.Forms.TabPage();
            this.tbTime = new System.Windows.Forms.TextBox();
            this.lblSeconds = new System.Windows.Forms.Label();
            this.tbMaxPressure = new System.Windows.Forms.TextBox();
            this.label1 = new System.Windows.Forms.Label();
            this.overpressureResultGrid = new System.Windows.Forms.DataGridView();
            this.outputWarning = new System.Windows.Forms.Label();
            this.notionalNozzleSelector = new System.Windows.Forms.ComboBox();
            this.label3 = new System.Windows.Forms.Label();
            this.colTime = new System.Windows.Forms.DataGridViewTextBoxColumn();
            this.colPressure = new System.Windows.Forms.DataGridViewTextBoxColumn();
            this.colDepth2 = new System.Windows.Forms.DataGridViewTextBoxColumn();
            this.colConcentration2 = new System.Windows.Forms.DataGridViewTextBoxColumn();
            this.pbPressure = new SandiaNationalLaboratories.Hyram.PictureBoxWithSave();
            this.pbFlammableMass = new SandiaNationalLaboratories.Hyram.PictureBoxWithSave();
            this.pbLayer = new SandiaNationalLaboratories.Hyram.PictureBoxWithSave();
            this.pbTrajectory = new SandiaNationalLaboratories.Hyram.PictureBoxWithSave();
            this.IOTabs.SuspendLayout();
            this.InputTab.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.OverpressureSpinner)).BeginInit();
            this.InputTabs.SuspendLayout();
            this.IndoorReleaseTab.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.GeometryPicture)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.InputGrid)).BeginInit();
            this.OutputOptionsTab.SuspendLayout();
            this.PlottingOptionsGroupBox.SuspendLayout();
            this.PressuresPerTimeGroupBox.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.PressuresPerTimeGrid)).BeginInit();
            this.PressureLinesGroupBox.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.PressureLinesGrid)).BeginInit();
            this.outputTab.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.splitContainer1)).BeginInit();
            this.splitContainer1.Panel1.SuspendLayout();
            this.splitContainer1.Panel2.SuspendLayout();
            this.splitContainer1.SuspendLayout();
            this.tcOutput.SuspendLayout();
            this.tpPressure.SuspendLayout();
            this.tpFlammableMass.SuspendLayout();
            this.tpLayer.SuspendLayout();
            this.tpTrajectory.SuspendLayout();
            this.tpData.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.overpressureResultGrid)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.pbPressure)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.pbFlammableMass)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.pbLayer)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.pbTrajectory)).BeginInit();
            this.SuspendLayout();
            // 
            // IOTabs
            // 
            this.IOTabs.Controls.Add(this.InputTab);
            this.IOTabs.Controls.Add(this.outputTab);
            this.IOTabs.Dock = System.Windows.Forms.DockStyle.Fill;
            this.IOTabs.Location = new System.Drawing.Point(0, 0);
            this.IOTabs.Name = "IOTabs";
            this.IOTabs.SelectedIndex = 0;
            this.IOTabs.Size = new System.Drawing.Size(938, 583);
            this.IOTabs.TabIndex = 0;
            // 
            // InputTab
            // 
            this.InputTab.Controls.Add(this.notionalNozzleSelector);
            this.InputTab.Controls.Add(this.label3);
            this.InputTab.Controls.Add(this.PhaseSelection);
            this.InputTab.Controls.Add(this.PhaseLabel);
            this.InputTab.Controls.Add(this.OverpressureSpinner);
            this.InputTab.Controls.Add(this.InputWarning);
            this.InputTab.Controls.Add(this.InputTabs);
            this.InputTab.Controls.Add(this.ExecuteBtn);
            this.InputTab.Location = new System.Drawing.Point(4, 22);
            this.InputTab.Name = "InputTab";
            this.InputTab.Padding = new System.Windows.Forms.Padding(3);
            this.InputTab.Size = new System.Drawing.Size(930, 557);
            this.InputTab.TabIndex = 0;
            this.InputTab.Text = "Input";
            this.InputTab.UseVisualStyleBackColor = true;
            // 
            // PhaseSelection
            // 
            this.PhaseSelection.DropDownStyle = System.Windows.Forms.ComboBoxStyle.DropDownList;
            this.PhaseSelection.DropDownWidth = 170;
            this.PhaseSelection.FormattingEnabled = true;
            this.PhaseSelection.Location = new System.Drawing.Point(167, 40);
            this.PhaseSelection.Name = "PhaseSelection";
            this.PhaseSelection.Size = new System.Drawing.Size(176, 21);
            this.PhaseSelection.TabIndex = 58;
            this.PhaseSelection.SelectionChangeCommitted += new System.EventHandler(this.PhaseSelection_SelectionChangeCommitted);
            // 
            // PhaseLabel
            // 
            this.PhaseLabel.AutoSize = true;
            this.PhaseLabel.Font = new System.Drawing.Font("Microsoft Sans Serif", 9F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.PhaseLabel.Location = new System.Drawing.Point(7, 41);
            this.PhaseLabel.Name = "PhaseLabel";
            this.PhaseLabel.Size = new System.Drawing.Size(71, 15);
            this.PhaseLabel.TabIndex = 57;
            this.PhaseLabel.Text = "Fluid phase";
            // 
            // OverpressureSpinner
            // 
            this.OverpressureSpinner.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Right)));
            this.OverpressureSpinner.Image = global::SandiaNationalLaboratories.Hyram.Properties.Resources.AjaxSpinner;
            this.OverpressureSpinner.Location = new System.Drawing.Point(810, 528);
            this.OverpressureSpinner.MinimumSize = new System.Drawing.Size(20, 20);
            this.OverpressureSpinner.Name = "OverpressureSpinner";
            this.OverpressureSpinner.Size = new System.Drawing.Size(32, 23);
            this.OverpressureSpinner.SizeMode = System.Windows.Forms.PictureBoxSizeMode.CenterImage;
            this.OverpressureSpinner.TabIndex = 18;
            this.OverpressureSpinner.TabStop = false;
            // 
            // InputWarning
            // 
            this.InputWarning.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Left)));
            this.InputWarning.AutoSize = true;
            this.InputWarning.BackColor = System.Drawing.Color.MistyRose;
            this.InputWarning.Font = new System.Drawing.Font("Microsoft Sans Serif", 9F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.InputWarning.ForeColor = System.Drawing.Color.Maroon;
            this.InputWarning.Location = new System.Drawing.Point(7, 527);
            this.InputWarning.Name = "InputWarning";
            this.InputWarning.Padding = new System.Windows.Forms.Padding(4);
            this.InputWarning.Size = new System.Drawing.Size(278, 23);
            this.InputWarning.TabIndex = 17;
            this.InputWarning.Text = "No errors noted and this label is invisible";
            this.InputWarning.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            this.InputWarning.Visible = false;
            // 
            // InputTabs
            // 
            this.InputTabs.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.InputTabs.Controls.Add(this.IndoorReleaseTab);
            this.InputTabs.Controls.Add(this.OutputOptionsTab);
            this.InputTabs.Location = new System.Drawing.Point(3, 78);
            this.InputTabs.Name = "InputTabs";
            this.InputTabs.SelectedIndex = 0;
            this.InputTabs.Size = new System.Drawing.Size(924, 444);
            this.InputTabs.TabIndex = 1;
            // 
            // IndoorReleaseTab
            // 
            this.IndoorReleaseTab.AutoScroll = true;
            this.IndoorReleaseTab.Controls.Add(this.GeometryPicture);
            this.IndoorReleaseTab.Controls.Add(this.InputGrid);
            this.IndoorReleaseTab.Location = new System.Drawing.Point(4, 22);
            this.IndoorReleaseTab.Name = "IndoorReleaseTab";
            this.IndoorReleaseTab.Padding = new System.Windows.Forms.Padding(3);
            this.IndoorReleaseTab.Size = new System.Drawing.Size(916, 418);
            this.IndoorReleaseTab.TabIndex = 0;
            this.IndoorReleaseTab.Text = "Indoor Release Parameters";
            this.IndoorReleaseTab.UseVisualStyleBackColor = true;
            // 
            // GeometryPicture
            // 
            this.GeometryPicture.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.GeometryPicture.BackgroundImage = ((System.Drawing.Image)(resources.GetObject("GeometryPicture.BackgroundImage")));
            this.GeometryPicture.BackgroundImageLayout = System.Windows.Forms.ImageLayout.Zoom;
            this.GeometryPicture.Location = new System.Drawing.Point(539, 6);
            this.GeometryPicture.Name = "GeometryPicture";
            this.GeometryPicture.Size = new System.Drawing.Size(255, 526);
            this.GeometryPicture.TabIndex = 1;
            this.GeometryPicture.TabStop = false;
            // 
            // InputGrid
            // 
            this.InputGrid.AllowUserToAddRows = false;
            this.InputGrid.AllowUserToDeleteRows = false;
            this.InputGrid.AutoSizeColumnsMode = System.Windows.Forms.DataGridViewAutoSizeColumnsMode.Fill;
            this.InputGrid.ColumnHeadersHeightSizeMode = System.Windows.Forms.DataGridViewColumnHeadersHeightSizeMode.AutoSize;
            this.InputGrid.Location = new System.Drawing.Point(3, 3);
            this.InputGrid.Name = "InputGrid";
            this.InputGrid.Size = new System.Drawing.Size(530, 484);
            this.InputGrid.TabIndex = 0;
            this.InputGrid.CellValueChanged += new System.Windows.Forms.DataGridViewCellEventHandler(this.InputGrid_CellValueChanged);
            // 
            // OutputOptionsTab
            // 
            this.OutputOptionsTab.Controls.Add(this.MaxTimeInput);
            this.OutputOptionsTab.Controls.Add(this.MaxTimeLabel);
            this.OutputOptionsTab.Controls.Add(this.PlottingOptionsGroupBox);
            this.OutputOptionsTab.Controls.Add(this.PlotTimesInput);
            this.OutputOptionsTab.Controls.Add(this.PlotTimesLabel);
            this.OutputOptionsTab.Location = new System.Drawing.Point(4, 22);
            this.OutputOptionsTab.Name = "OutputOptionsTab";
            this.OutputOptionsTab.Padding = new System.Windows.Forms.Padding(3);
            this.OutputOptionsTab.Size = new System.Drawing.Size(916, 418);
            this.OutputOptionsTab.TabIndex = 1;
            this.OutputOptionsTab.Text = "Output Options";
            this.OutputOptionsTab.UseVisualStyleBackColor = true;
            // 
            // MaxTimeInput
            // 
            this.MaxTimeInput.Location = new System.Drawing.Point(230, 39);
            this.MaxTimeInput.Name = "MaxTimeInput";
            this.MaxTimeInput.Size = new System.Drawing.Size(100, 20);
            this.MaxTimeInput.TabIndex = 6;
            this.MaxTimeInput.TextChanged += new System.EventHandler(this.MaxTimeInput_TextChanged);
            // 
            // MaxTimeLabel
            // 
            this.MaxTimeLabel.AutoSize = true;
            this.MaxTimeLabel.Location = new System.Drawing.Point(6, 42);
            this.MaxTimeLabel.Name = "MaxTimeLabel";
            this.MaxTimeLabel.Size = new System.Drawing.Size(136, 13);
            this.MaxTimeLabel.TabIndex = 5;
            this.MaxTimeLabel.Text = "Maximum time (in seconds):";
            // 
            // PlottingOptionsGroupBox
            // 
            this.PlottingOptionsGroupBox.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.PlottingOptionsGroupBox.Controls.Add(this.PressuresPerTimeCheckbox);
            this.PlottingOptionsGroupBox.Controls.Add(this.PressuresPerTimeGroupBox);
            this.PlottingOptionsGroupBox.Controls.Add(this.PressureLinesGroupBox);
            this.PlottingOptionsGroupBox.Controls.Add(this.PressureLinesCheckbox);
            this.PlottingOptionsGroupBox.Location = new System.Drawing.Point(6, 66);
            this.PlottingOptionsGroupBox.Name = "PlottingOptionsGroupBox";
            this.PlottingOptionsGroupBox.Size = new System.Drawing.Size(904, 349);
            this.PlottingOptionsGroupBox.TabIndex = 4;
            this.PlottingOptionsGroupBox.TabStop = false;
            this.PlottingOptionsGroupBox.Text = "Plotting options";
            // 
            // PressuresPerTimeCheckbox
            // 
            this.PressuresPerTimeCheckbox.AutoSize = true;
            this.PressuresPerTimeCheckbox.Checked = true;
            this.PressuresPerTimeCheckbox.CheckState = System.Windows.Forms.CheckState.Checked;
            this.PressuresPerTimeCheckbox.Cursor = System.Windows.Forms.Cursors.Arrow;
            this.PressuresPerTimeCheckbox.Location = new System.Drawing.Point(287, 16);
            this.PressuresPerTimeCheckbox.Name = "PressuresPerTimeCheckbox";
            this.PressuresPerTimeCheckbox.Size = new System.Drawing.Size(186, 17);
            this.PressuresPerTimeCheckbox.TabIndex = 10;
            this.PressuresPerTimeCheckbox.Text = "Mark chart with pressures at times";
            this.PressuresPerTimeCheckbox.UseVisualStyleBackColor = true;
            this.PressuresPerTimeCheckbox.CheckedChanged += new System.EventHandler(this.PressuresPerTimeCheckbox_CheckedChanged);
            // 
            // PressuresPerTimeGroupBox
            // 
            this.PressuresPerTimeGroupBox.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left)));
            this.PressuresPerTimeGroupBox.Controls.Add(this.PressuresPerTimeGrid);
            this.PressuresPerTimeGroupBox.Location = new System.Drawing.Point(284, 40);
            this.PressuresPerTimeGroupBox.Name = "PressuresPerTimeGroupBox";
            this.PressuresPerTimeGroupBox.Size = new System.Drawing.Size(316, 300);
            this.PressuresPerTimeGroupBox.TabIndex = 9;
            this.PressuresPerTimeGroupBox.TabStop = false;
            this.PressuresPerTimeGroupBox.Text = "Place dots where pressure/time intersect";
            // 
            // PressuresPerTimeGrid
            // 
            this.PressuresPerTimeGrid.ColumnHeadersHeightSizeMode = System.Windows.Forms.DataGridViewColumnHeadersHeightSizeMode.AutoSize;
            this.PressuresPerTimeGrid.Columns.AddRange(new System.Windows.Forms.DataGridViewColumn[] {
            this.MarkTimeGridCol,
            this.MarkPressureGridCol});
            this.PressuresPerTimeGrid.Dock = System.Windows.Forms.DockStyle.Fill;
            this.PressuresPerTimeGrid.Location = new System.Drawing.Point(3, 16);
            this.PressuresPerTimeGrid.Name = "PressuresPerTimeGrid";
            this.PressuresPerTimeGrid.Size = new System.Drawing.Size(310, 281);
            this.PressuresPerTimeGrid.TabIndex = 4;
            this.PressuresPerTimeGrid.CellValueChanged += new System.Windows.Forms.DataGridViewCellEventHandler(this.PressuresPerTimeGrid_CellValueChanged);
            this.PressuresPerTimeGrid.RowsRemoved += new System.Windows.Forms.DataGridViewRowsRemovedEventHandler(this.PressuresPerTimeGrid_RowsRemoved);
            this.PressuresPerTimeGrid.SortCompare += new System.Windows.Forms.DataGridViewSortCompareEventHandler(this.PressuresPerTimeGrid_SortCompare);
            // 
            // MarkTimeGridCol
            // 
            this.MarkTimeGridCol.HeaderText = "Time (Seconds)";
            this.MarkTimeGridCol.Name = "MarkTimeGridCol";
            // 
            // MarkPressureGridCol
            // 
            this.MarkPressureGridCol.HeaderText = "Pressure (kPa)";
            this.MarkPressureGridCol.Name = "MarkPressureGridCol";
            // 
            // PressureLinesGroupBox
            // 
            this.PressureLinesGroupBox.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left)));
            this.PressureLinesGroupBox.Controls.Add(this.PressureLinesGrid);
            this.PressureLinesGroupBox.Location = new System.Drawing.Point(9, 37);
            this.PressureLinesGroupBox.Name = "PressureLinesGroupBox";
            this.PressureLinesGroupBox.Size = new System.Drawing.Size(225, 306);
            this.PressureLinesGroupBox.TabIndex = 8;
            this.PressureLinesGroupBox.TabStop = false;
            this.PressureLinesGroupBox.Text = "Specify pressures in kPa";
            // 
            // PressureLinesGrid
            // 
            this.PressureLinesGrid.ColumnHeadersHeightSizeMode = System.Windows.Forms.DataGridViewColumnHeadersHeightSizeMode.AutoSize;
            this.PressureLinesGrid.Columns.AddRange(new System.Windows.Forms.DataGridViewColumn[] {
            this.PressureLineCol});
            this.PressureLinesGrid.Dock = System.Windows.Forms.DockStyle.Fill;
            this.PressureLinesGrid.Location = new System.Drawing.Point(3, 16);
            this.PressureLinesGrid.Name = "PressureLinesGrid";
            this.PressureLinesGrid.Size = new System.Drawing.Size(219, 287);
            this.PressureLinesGrid.TabIndex = 6;
            this.PressureLinesGrid.CellValueChanged += new System.Windows.Forms.DataGridViewCellEventHandler(this.PressureLinesGrid_CellValueChanged);
            this.PressureLinesGrid.RowsRemoved += new System.Windows.Forms.DataGridViewRowsRemovedEventHandler(this.PressureLinesGrid_RowsRemoved);
            this.PressureLinesGrid.SortCompare += new System.Windows.Forms.DataGridViewSortCompareEventHandler(this.PressureLinesGrid_SortCompare);
            // 
            // PressureLineCol
            // 
            this.PressureLineCol.HeaderText = "Pressure (kPa)";
            this.PressureLineCol.Name = "PressureLineCol";
            // 
            // PressureLinesCheckbox
            // 
            this.PressureLinesCheckbox.AutoSize = true;
            this.PressureLinesCheckbox.Checked = true;
            this.PressureLinesCheckbox.CheckState = System.Windows.Forms.CheckState.Checked;
            this.PressureLinesCheckbox.Cursor = System.Windows.Forms.Cursors.Arrow;
            this.PressureLinesCheckbox.Location = new System.Drawing.Point(10, 16);
            this.PressureLinesCheckbox.Name = "PressureLinesCheckbox";
            this.PressureLinesCheckbox.Size = new System.Drawing.Size(270, 17);
            this.PressureLinesCheckbox.TabIndex = 7;
            this.PressureLinesCheckbox.Text = "Draw horizontal lines on chart at specified pressures";
            this.PressureLinesCheckbox.UseVisualStyleBackColor = true;
            this.PressureLinesCheckbox.CheckedChanged += new System.EventHandler(this.PressureLinesCheckbox_CheckedChanged);
            // 
            // PlotTimesInput
            // 
            this.PlotTimesInput.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.PlotTimesInput.Location = new System.Drawing.Point(230, 10);
            this.PlotTimesInput.Name = "PlotTimesInput";
            this.PlotTimesInput.Size = new System.Drawing.Size(668, 20);
            this.PlotTimesInput.TabIndex = 1;
            this.PlotTimesInput.TextChanged += new System.EventHandler(this.PlotTimesInput_TextChanged);
            // 
            // PlotTimesLabel
            // 
            this.PlotTimesLabel.AutoSize = true;
            this.PlotTimesLabel.Location = new System.Drawing.Point(6, 12);
            this.PlotTimesLabel.Name = "PlotTimesLabel";
            this.PlotTimesLabel.Size = new System.Drawing.Size(218, 13);
            this.PlotTimesLabel.TabIndex = 0;
            this.PlotTimesLabel.Text = "Output pressures at these times (in seconds):";
            // 
            // ExecuteBtn
            // 
            this.ExecuteBtn.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Right)));
            this.ExecuteBtn.Location = new System.Drawing.Point(848, 528);
            this.ExecuteBtn.Name = "ExecuteBtn";
            this.ExecuteBtn.Size = new System.Drawing.Size(75, 23);
            this.ExecuteBtn.TabIndex = 16;
            this.ExecuteBtn.Text = "Calculate";
            this.ExecuteBtn.UseVisualStyleBackColor = true;
            this.ExecuteBtn.Click += new System.EventHandler(this.ExecuteBtn_Click);
            // 
            // outputTab
            // 
            this.outputTab.Controls.Add(this.splitContainer1);
            this.outputTab.Location = new System.Drawing.Point(4, 22);
            this.outputTab.Name = "outputTab";
            this.outputTab.Padding = new System.Windows.Forms.Padding(3);
            this.outputTab.Size = new System.Drawing.Size(930, 557);
            this.outputTab.TabIndex = 1;
            this.outputTab.Text = "Output";
            this.outputTab.UseVisualStyleBackColor = true;
            // 
            // splitContainer1
            // 
            this.splitContainer1.Dock = System.Windows.Forms.DockStyle.Fill;
            this.splitContainer1.Location = new System.Drawing.Point(3, 3);
            this.splitContainer1.Name = "splitContainer1";
            this.splitContainer1.Orientation = System.Windows.Forms.Orientation.Horizontal;
            // 
            // splitContainer1.Panel1
            // 
            this.splitContainer1.Panel1.Controls.Add(this.tcOutput);
            // 
            // splitContainer1.Panel2
            // 
            this.splitContainer1.Panel2.Controls.Add(this.outputWarning);
            this.splitContainer1.Size = new System.Drawing.Size(924, 551);
            this.splitContainer1.SplitterDistance = 522;
            this.splitContainer1.TabIndex = 19;
            // 
            // tcOutput
            // 
            this.tcOutput.Controls.Add(this.tpPressure);
            this.tcOutput.Controls.Add(this.tpFlammableMass);
            this.tcOutput.Controls.Add(this.tpLayer);
            this.tcOutput.Controls.Add(this.tpTrajectory);
            this.tcOutput.Controls.Add(this.tpData);
            this.tcOutput.Dock = System.Windows.Forms.DockStyle.Fill;
            this.tcOutput.Location = new System.Drawing.Point(0, 0);
            this.tcOutput.Name = "tcOutput";
            this.tcOutput.SelectedIndex = 0;
            this.tcOutput.Size = new System.Drawing.Size(924, 522);
            this.tcOutput.TabIndex = 3;
            // 
            // tpPressure
            // 
            this.tpPressure.Controls.Add(this.label2);
            this.tpPressure.Controls.Add(this.pbPressure);
            this.tpPressure.Location = new System.Drawing.Point(4, 22);
            this.tpPressure.Name = "tpPressure";
            this.tpPressure.Padding = new System.Windows.Forms.Padding(3);
            this.tpPressure.Size = new System.Drawing.Size(916, 496);
            this.tpPressure.TabIndex = 0;
            this.tpPressure.Text = "Pressure plot";
            this.tpPressure.UseVisualStyleBackColor = true;
            // 
            // label2
            // 
            this.label2.Dock = System.Windows.Forms.DockStyle.Top;
            this.label2.Location = new System.Drawing.Point(3, 3);
            this.label2.Name = "label2";
            this.label2.Padding = new System.Windows.Forms.Padding(3);
            this.label2.Size = new System.Drawing.Size(910, 24);
            this.label2.TabIndex = 1;
            this.label2.Text = "Time-history of the overpressure that would develop if the accumulated hydrogen w" +
    "ere to be ignited after some delay after the leak started";
            // 
            // tpFlammableMass
            // 
            this.tpFlammableMass.Controls.Add(this.pbFlammableMass);
            this.tpFlammableMass.Location = new System.Drawing.Point(4, 22);
            this.tpFlammableMass.Name = "tpFlammableMass";
            this.tpFlammableMass.Size = new System.Drawing.Size(916, 496);
            this.tpFlammableMass.TabIndex = 4;
            this.tpFlammableMass.Text = "Flammable mass plot";
            this.tpFlammableMass.UseVisualStyleBackColor = true;
            // 
            // tpLayer
            // 
            this.tpLayer.Controls.Add(this.pbLayer);
            this.tpLayer.Location = new System.Drawing.Point(4, 22);
            this.tpLayer.Name = "tpLayer";
            this.tpLayer.Padding = new System.Windows.Forms.Padding(3);
            this.tpLayer.Size = new System.Drawing.Size(916, 496);
            this.tpLayer.TabIndex = 1;
            this.tpLayer.Text = "Layer plot";
            this.tpLayer.UseVisualStyleBackColor = true;
            // 
            // tpTrajectory
            // 
            this.tpTrajectory.Controls.Add(this.pbTrajectory);
            this.tpTrajectory.Location = new System.Drawing.Point(4, 22);
            this.tpTrajectory.Name = "tpTrajectory";
            this.tpTrajectory.Size = new System.Drawing.Size(916, 496);
            this.tpTrajectory.TabIndex = 3;
            this.tpTrajectory.Text = "Trajectory plot";
            this.tpTrajectory.UseVisualStyleBackColor = true;
            // 
            // tpData
            // 
            this.tpData.Controls.Add(this.tbTime);
            this.tpData.Controls.Add(this.lblSeconds);
            this.tpData.Controls.Add(this.tbMaxPressure);
            this.tpData.Controls.Add(this.label1);
            this.tpData.Controls.Add(this.overpressureResultGrid);
            this.tpData.Location = new System.Drawing.Point(4, 22);
            this.tpData.Name = "tpData";
            this.tpData.Size = new System.Drawing.Size(916, 496);
            this.tpData.TabIndex = 2;
            this.tpData.Text = "Data";
            this.tpData.UseVisualStyleBackColor = true;
            // 
            // tbTime
            // 
            this.tbTime.Enabled = false;
            this.tbTime.Location = new System.Drawing.Point(172, 38);
            this.tbTime.Name = "tbTime";
            this.tbTime.Size = new System.Drawing.Size(184, 20);
            this.tbTime.TabIndex = 6;
            // 
            // lblSeconds
            // 
            this.lblSeconds.AutoSize = true;
            this.lblSeconds.Location = new System.Drawing.Point(16, 41);
            this.lblSeconds.Name = "lblSeconds";
            this.lblSeconds.Size = new System.Drawing.Size(150, 13);
            this.lblSeconds.TabIndex = 5;
            this.lblSeconds.Text = "Time this Occurred (Seconds):";
            // 
            // tbMaxPressure
            // 
            this.tbMaxPressure.Enabled = false;
            this.tbMaxPressure.Location = new System.Drawing.Point(172, 10);
            this.tbMaxPressure.Name = "tbMaxPressure";
            this.tbMaxPressure.Size = new System.Drawing.Size(184, 20);
            this.tbMaxPressure.TabIndex = 4;
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(16, 17);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(142, 13);
            this.label1.TabIndex = 3;
            this.label1.Text = "Maximum Overpressure (Pa):";
            // 
            // overpressureResultGrid
            // 
            this.overpressureResultGrid.AllowUserToAddRows = false;
            this.overpressureResultGrid.AllowUserToDeleteRows = false;
            this.overpressureResultGrid.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left)));
            this.overpressureResultGrid.AutoSizeColumnsMode = System.Windows.Forms.DataGridViewAutoSizeColumnsMode.Fill;
            dataGridViewCellStyle2.Alignment = System.Windows.Forms.DataGridViewContentAlignment.MiddleLeft;
            dataGridViewCellStyle2.BackColor = System.Drawing.SystemColors.Control;
            dataGridViewCellStyle2.Font = new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            dataGridViewCellStyle2.ForeColor = System.Drawing.SystemColors.WindowText;
            dataGridViewCellStyle2.SelectionBackColor = System.Drawing.SystemColors.Highlight;
            dataGridViewCellStyle2.SelectionForeColor = System.Drawing.SystemColors.HighlightText;
            dataGridViewCellStyle2.WrapMode = System.Windows.Forms.DataGridViewTriState.True;
            this.overpressureResultGrid.ColumnHeadersDefaultCellStyle = dataGridViewCellStyle2;
            this.overpressureResultGrid.ColumnHeadersHeightSizeMode = System.Windows.Forms.DataGridViewColumnHeadersHeightSizeMode.AutoSize;
            this.overpressureResultGrid.Columns.AddRange(new System.Windows.Forms.DataGridViewColumn[] {
            this.colTime,
            this.colPressure,
            this.colDepth2,
            this.colConcentration2});
            this.overpressureResultGrid.Location = new System.Drawing.Point(19, 73);
            this.overpressureResultGrid.Name = "overpressureResultGrid";
            this.overpressureResultGrid.RowHeadersVisible = false;
            this.overpressureResultGrid.Size = new System.Drawing.Size(536, 399);
            this.overpressureResultGrid.TabIndex = 2;
            // 
            // outputWarning
            // 
            this.outputWarning.BackColor = System.Drawing.Color.PaleGoldenrod;
            this.outputWarning.Dock = System.Windows.Forms.DockStyle.Fill;
            this.outputWarning.Font = new System.Drawing.Font("Microsoft Sans Serif", 9F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.outputWarning.ForeColor = System.Drawing.Color.DarkGoldenrod;
            this.outputWarning.Location = new System.Drawing.Point(0, 0);
            this.outputWarning.Name = "outputWarning";
            this.outputWarning.Size = new System.Drawing.Size(924, 25);
            this.outputWarning.TabIndex = 18;
            this.outputWarning.Text = "blank";
            this.outputWarning.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            this.outputWarning.Visible = false;
            // 
            // notionalNozzleSelector
            // 
            this.notionalNozzleSelector.DropDownStyle = System.Windows.Forms.ComboBoxStyle.DropDownList;
            this.notionalNozzleSelector.DropDownWidth = 170;
            this.notionalNozzleSelector.FormattingEnabled = true;
            this.notionalNozzleSelector.Location = new System.Drawing.Point(167, 13);
            this.notionalNozzleSelector.Name = "notionalNozzleSelector";
            this.notionalNozzleSelector.Size = new System.Drawing.Size(176, 21);
            this.notionalNozzleSelector.TabIndex = 60;
            this.notionalNozzleSelector.SelectionChangeCommitted += new System.EventHandler(this.notionalNozzleSelector_SelectionChangeCommitted);
            // 
            // label3
            // 
            this.label3.Font = new System.Drawing.Font("Microsoft Sans Serif", 9F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.label3.Location = new System.Drawing.Point(7, 14);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(153, 15);
            this.label3.TabIndex = 59;
            this.label3.Text = "Notional nozzle model:";
            // 
            // colTime
            // 
            this.colTime.HeaderText = "Time (Seconds)";
            this.colTime.Name = "colTime";
            this.colTime.ReadOnly = true;
            // 
            // colPressure
            // 
            this.colPressure.HeaderText = "Combined Pressure (Pa)";
            this.colPressure.Name = "colPressure";
            this.colPressure.ReadOnly = true;
            // 
            // colDepth2
            // 
            this.colDepth2.HeaderText = "Depth (m)";
            this.colDepth2.Name = "colDepth2";
            this.colDepth2.ReadOnly = true;
            // 
            // colConcentration2
            // 
            this.colConcentration2.HeaderText = "Concentration (%)";
            this.colConcentration2.Name = "colConcentration2";
            this.colConcentration2.ReadOnly = true;
            // 
            // pbPressure
            // 
            this.pbPressure.Dock = System.Windows.Forms.DockStyle.Bottom;
            this.pbPressure.Location = new System.Drawing.Point(3, 30);
            this.pbPressure.Name = "pbPressure";
            this.pbPressure.Size = new System.Drawing.Size(910, 463);
            this.pbPressure.TabIndex = 0;
            this.pbPressure.TabStop = false;
            // 
            // pbFlammableMass
            // 
            this.pbFlammableMass.Dock = System.Windows.Forms.DockStyle.Fill;
            this.pbFlammableMass.Location = new System.Drawing.Point(0, 0);
            this.pbFlammableMass.Name = "pbFlammableMass";
            this.pbFlammableMass.Size = new System.Drawing.Size(916, 496);
            this.pbFlammableMass.TabIndex = 1;
            this.pbFlammableMass.TabStop = false;
            // 
            // pbLayer
            // 
            this.pbLayer.Dock = System.Windows.Forms.DockStyle.Fill;
            this.pbLayer.Location = new System.Drawing.Point(3, 3);
            this.pbLayer.Name = "pbLayer";
            this.pbLayer.Size = new System.Drawing.Size(910, 490);
            this.pbLayer.TabIndex = 0;
            this.pbLayer.TabStop = false;
            // 
            // pbTrajectory
            // 
            this.pbTrajectory.Dock = System.Windows.Forms.DockStyle.Fill;
            this.pbTrajectory.Location = new System.Drawing.Point(0, 0);
            this.pbTrajectory.Name = "pbTrajectory";
            this.pbTrajectory.Size = new System.Drawing.Size(916, 496);
            this.pbTrajectory.TabIndex = 1;
            this.pbTrajectory.TabStop = false;
            // 
            // IndoorReleaseForm
            // 
            this.Controls.Add(this.IOTabs);
            this.Name = "IndoorReleaseForm";
            this.Size = new System.Drawing.Size(938, 583);
            this.Load += new System.EventHandler(this.IndoorReleaseForm_Load);
            this.IOTabs.ResumeLayout(false);
            this.InputTab.ResumeLayout(false);
            this.InputTab.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.OverpressureSpinner)).EndInit();
            this.InputTabs.ResumeLayout(false);
            this.IndoorReleaseTab.ResumeLayout(false);
            ((System.ComponentModel.ISupportInitialize)(this.GeometryPicture)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.InputGrid)).EndInit();
            this.OutputOptionsTab.ResumeLayout(false);
            this.OutputOptionsTab.PerformLayout();
            this.PlottingOptionsGroupBox.ResumeLayout(false);
            this.PlottingOptionsGroupBox.PerformLayout();
            this.PressuresPerTimeGroupBox.ResumeLayout(false);
            ((System.ComponentModel.ISupportInitialize)(this.PressuresPerTimeGrid)).EndInit();
            this.PressureLinesGroupBox.ResumeLayout(false);
            ((System.ComponentModel.ISupportInitialize)(this.PressureLinesGrid)).EndInit();
            this.outputTab.ResumeLayout(false);
            this.splitContainer1.Panel1.ResumeLayout(false);
            this.splitContainer1.Panel2.ResumeLayout(false);
            ((System.ComponentModel.ISupportInitialize)(this.splitContainer1)).EndInit();
            this.splitContainer1.ResumeLayout(false);
            this.tcOutput.ResumeLayout(false);
            this.tpPressure.ResumeLayout(false);
            this.tpFlammableMass.ResumeLayout(false);
            this.tpLayer.ResumeLayout(false);
            this.tpTrajectory.ResumeLayout(false);
            this.tpData.ResumeLayout(false);
            this.tpData.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.overpressureResultGrid)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.pbPressure)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.pbFlammableMass)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.pbLayer)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.pbTrajectory)).EndInit();
            this.ResumeLayout(false);

		}

		#endregion

		private TabControl IOTabs;
		private TabPage InputTab;
		private DataGridView InputGrid;
		private TabPage outputTab;
		private Button ExecuteBtn;
		private TabControl InputTabs;
		private TabPage IndoorReleaseTab;
		private TabPage OutputOptionsTab;
		private TextBox PlotTimesInput;
		private Label PlotTimesLabel;
		private GroupBox PlottingOptionsGroupBox;
		private DataGridView PressuresPerTimeGrid;
		private DataGridView PressureLinesGrid;
		private DataGridViewTextBoxColumn PressureLineCol;
		private DataGridViewTextBoxColumn MarkTimeGridCol;
		private DataGridViewTextBoxColumn MarkPressureGridCol;
		private DataGridView overpressureResultGrid;
		private CheckBox PressureLinesCheckbox;
		private GroupBox PressureLinesGroupBox;
		private GroupBox PressuresPerTimeGroupBox;
		private CheckBox PressuresPerTimeCheckbox;
		private PictureBox GeometryPicture;
		private TabControl tcOutput;
		private TabPage tpPressure;
		private PictureBoxWithSave pbPressure;
		private TabPage tpLayer;
		private PictureBoxWithSave pbLayer;
		private TabPage tpData;
		private TextBox tbTime;
		private Label lblSeconds;
		private TextBox tbMaxPressure;
		private Label label1;
		private Label InputWarning;
		private TextBox MaxTimeInput;
		private Label MaxTimeLabel;
		private TabPage tpTrajectory;
		private PictureBoxWithSave pbTrajectory;
		private TabPage tpFlammableMass;
		private PictureBoxWithSave pbFlammableMass;
        private PictureBox OverpressureSpinner;
        private Label outputWarning;
        private SplitContainer splitContainer1;
        private ComboBox PhaseSelection;
        private Label PhaseLabel;
        private Label label2;
        private ComboBox notionalNozzleSelector;
        private Label label3;
        private DataGridViewTextBoxColumn colTime;
        private DataGridViewTextBoxColumn colPressure;
        private DataGridViewTextBoxColumn colDepth2;
        private DataGridViewTextBoxColumn colConcentration2;
    }
}
