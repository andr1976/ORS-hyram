/*
Copyright 2015-2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, the U.S.Government retains certain
rights in this software.

You should have received a copy of the GNU General Public License along with
HyRAM+. If not, see https://www.gnu.org/licenses/.
*/

using System;
using System.ComponentModel;
using System.Diagnostics;
using System.IO;
using System.Windows.Forms;
using SandiaNationalLaboratories.Hyram.Resources;

namespace SandiaNationalLaboratories.Hyram
{
    public class AboutForm : Form
    {
        private readonly Container _components = null;

        private Button _okBtn;

        private Label _buildDayText;
        private Label _buildByText;
        private Label _versionText;
        private LinkLabel _weblink;

        private readonly string _mWebsiteUrl;
        private PictureBox _bannerLogo;
        private RichTextBox _copyrightText;

        public AboutForm()
        {
            InitializeComponent();

            _mWebsiteUrl = "http://hyram.sandia.gov";
            //Text = "About " + Application.ProductName;
            _versionText.Text = "Version " + CreateVersionString();

            var buildDateTime = File.GetLastWriteTime(Application.ExecutablePath);
            _buildDayText.Text = "Built on " + buildDateTime.ToLongDateString() + " at " +
                               buildDateTime.ToLongTimeString();

            _copyrightText.Rtf = Narratives.CopyrightFormatted;
        }

        private string CreateVersionString()
        {
            return Application.ProductVersion;
            //var components = Application.ProductVersion.Split('.');
            //string result = null;

            //for (var index = 0; index < components.Length; index++)
            //    if (index == components.Length - 1)
            //    {
            //        //result += ", build " + components[index];
            //    }
            //    else
            //    {
            //        if (result == null)
            //            result = components[index];
            //        else
            //            result += "." + components[index];
            //    }

            //return result;
        }

        /// <summary>
        ///     Clean up any resources being used.
        /// </summary>
        protected override void Dispose(bool disposing)
        {
            if (disposing)
                if (_components != null)
                    _components.Dispose();

            base.Dispose(disposing);
        }

#region Windows Form Designer generated code

        /// <summary>
        ///     Required method for Designer support - do not modify
        ///     the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            System.ComponentModel.ComponentResourceManager resources = new System.ComponentModel.ComponentResourceManager(typeof(AboutForm));
            this._versionText = new System.Windows.Forms.Label();
            this._buildByText = new System.Windows.Forms.Label();
            this._okBtn = new System.Windows.Forms.Button();
            this._buildDayText = new System.Windows.Forms.Label();
            this._weblink = new System.Windows.Forms.LinkLabel();
            this._copyrightText = new System.Windows.Forms.RichTextBox();
            this._bannerLogo = new System.Windows.Forms.PictureBox();
            ((System.ComponentModel.ISupportInitialize)(this._bannerLogo)).BeginInit();
            this.SuspendLayout();
            // 
            // _versionText
            // 
            this._versionText.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this._versionText.Location = new System.Drawing.Point(15, 93);
            this._versionText.Name = "_versionText";
            this._versionText.Size = new System.Drawing.Size(376, 26);
            this._versionText.TabIndex = 1;
            this._versionText.Text = "Version XX.XXX.XXX";
            this._versionText.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            // 
            // _buildByText
            // 
            this._buildByText.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this._buildByText.Location = new System.Drawing.Point(15, 129);
            this._buildByText.Name = "_buildByText";
            this._buildByText.Size = new System.Drawing.Size(376, 106);
            this._buildByText.TabIndex = 2;
            this._buildByText.Text = resources.GetString("_buildByText.Text");
            this._buildByText.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            // 
            // _okBtn
            // 
            this._okBtn.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Right)));
            this._okBtn.DialogResult = System.Windows.Forms.DialogResult.Cancel;
            this._okBtn.Location = new System.Drawing.Point(333, 371);
            this._okBtn.Name = "_okBtn";
            this._okBtn.Size = new System.Drawing.Size(67, 21);
            this._okBtn.TabIndex = 5;
            this._okBtn.Text = "OK";
            this._okBtn.Click += new System.EventHandler(this.OkBtn_Click);
            // 
            // _buildDayText
            // 
            this._buildDayText.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this._buildDayText.Location = new System.Drawing.Point(15, 113);
            this._buildDayText.Name = "_buildDayText";
            this._buildDayText.Size = new System.Drawing.Size(376, 16);
            this._buildDayText.TabIndex = 7;
            this._buildDayText.Text = "Built on";
            this._buildDayText.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            // 
            // _weblink
            // 
            this._weblink.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Left)));
            this._weblink.AutoSize = true;
            this._weblink.Location = new System.Drawing.Point(12, 375);
            this._weblink.Name = "_weblink";
            this._weblink.Size = new System.Drawing.Size(92, 13);
            this._weblink.TabIndex = 8;
            this._weblink.TabStop = true;
            this._weblink.Text = "HyRAM+ Website";
            this._weblink.LinkClicked += new System.Windows.Forms.LinkLabelLinkClickedEventHandler(this.Weblink_LinkClicked);
            // 
            // _copyrightText
            // 
            this._copyrightText.BackColor = System.Drawing.SystemColors.Control;
            this._copyrightText.BorderStyle = System.Windows.Forms.BorderStyle.None;
            this._copyrightText.Location = new System.Drawing.Point(18, 238);
            this._copyrightText.Name = "_copyrightText";
            this._copyrightText.ReadOnly = true;
            this._copyrightText.ScrollBars = System.Windows.Forms.RichTextBoxScrollBars.Vertical;
            this._copyrightText.Size = new System.Drawing.Size(378, 117);
            this._copyrightText.TabIndex = 11;
            this._copyrightText.Text = "";
            // 
            // _bannerLogo
            // 
            this._bannerLogo.Image = global::SandiaNationalLaboratories.Hyram.Properties.Resources.BannerLogov2Thin2;
            this._bannerLogo.Location = new System.Drawing.Point(0, -6);
            this._bannerLogo.Name = "_bannerLogo";
            this._bannerLogo.Size = new System.Drawing.Size(411, 105);
            this._bannerLogo.SizeMode = System.Windows.Forms.PictureBoxSizeMode.Zoom;
            this._bannerLogo.TabIndex = 9;
            this._bannerLogo.TabStop = false;
            // 
            // AboutForm
            // 
            this.AcceptButton = this._okBtn;
            this.AutoScaleBaseSize = new System.Drawing.Size(5, 13);
            this.BackgroundImageLayout = System.Windows.Forms.ImageLayout.Zoom;
            this.CancelButton = this._okBtn;
            this.ClientSize = new System.Drawing.Size(412, 399);
            this.Controls.Add(this._copyrightText);
            this.Controls.Add(this._bannerLogo);
            this.Controls.Add(this._weblink);
            this.Controls.Add(this._buildDayText);
            this.Controls.Add(this._okBtn);
            this.Controls.Add(this._buildByText);
            this.Controls.Add(this._versionText);
            this.FormBorderStyle = System.Windows.Forms.FormBorderStyle.FixedDialog;
            this.Icon = ((System.Drawing.Icon)(resources.GetObject("$this.Icon")));
            this.MaximizeBox = false;
            this.MinimizeBox = false;
            this.Name = "AboutForm";
            this.StartPosition = System.Windows.Forms.FormStartPosition.CenterScreen;
            this.Text = "About HyRAM+";
            ((System.ComponentModel.ISupportInitialize)(this._bannerLogo)).EndInit();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

#endregion

        private void OkBtn_Click(object sender, EventArgs e)
        {
            Close();
        }

        private void Weblink_LinkClicked(object sender, LinkLabelLinkClickedEventArgs e)
        {
            try
            {
                Process.Start(_mWebsiteUrl);
            }
            catch (Exception ex)
            {
                MessageBox.Show("Could not access website due to error: " + ex);
            }
        }
    }
}