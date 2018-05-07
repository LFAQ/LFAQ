using Microsoft.Win32;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Collections.ObjectModel;
using System.Collections;
using System.Diagnostics;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Windows.Forms.Design;
using System.Windows.Forms;
using System.Windows.Threading;



namespace LFAQ
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    /// 
    public partial class MainWindow : Window
    {
        public bool bIfQuanted;
        public bool bIfRunLoad;
        public bool bIfRunQuant;
        private bool bIfExportResult; //Make sure one run corresponds to one export;
        public bool bIfSaveStandProteins;  // Make sure the loading run after saving.
        private static System.Diagnostics.Process p;
        public bool bIfInstallR;
        private bool bIfLoadRunning;
        private bool bIfQuantiRunning;      
        Process LoadProcess=new Process();
        Process QuantproBach= new Process();
        public ObservableCollection<ProteinResult> ProteinsResult { get; set; }
        public ObservableCollection<StandProteinResult> StandProteinsResult { get; set; }
        public ObservableCollection<ExperimentName> SelectExperiments { get; set; }

        private delegate void ChangeStatusBar(String arg);
        private delegate void DelegateShowError();
        private delegate void ChangeParamIsEnable();
        public MainWindow()
        {           
            InitializeComponent();
            bIfQuanted = false;
            bIfExportResult = true; 
            bIfRunLoad = false;
            bIfRunQuant = false;
            bIfSaveStandProteins = false;
            bIfQuantiRunning = false;
            bIfLoadRunning = false;
            bIfInstallR = true;
            string strDefaultParameters = "Default.params";
            LoadParameters(strDefaultParameters);
            this.Closing += ClosingApplication;
        }

        private void ClosingApplication(object sender, CancelEventArgs e)
        {
            try
            {
                MessageBoxResult result = System.Windows.MessageBox.Show("Are you sure to exit？", "Closing", MessageBoxButton.YesNo, MessageBoxImage.Question);
                //Close the windown
                if (result == MessageBoxResult.Yes)
                {
                    KillProcesses();
                    e.Cancel = false;
                }
                //cancel the decision to close the window
                if (result == MessageBoxResult.No)
                {
                    e.Cancel = true;
                }
            }
            catch (NotImplementedException nie)
            {
                //throw new NotImplementedException();
            }
            catch(Exception ex)
            {
                System.Windows.MessageBox.Show(ex.ToString());              
            }

        }
        #region Menu Handlers
        protected void FileExit_Click(object sender, RoutedEventArgs args)
        {
            // Close this window.
            KillProcesses();
            this.Close();
        }
        private void UpdateStatusBar(String args)
        {
            statBarText.Text = args;
        }
        
        //protected void MouseEnterExitArea(object sender, RoutedEventArgs args)
        //{
        //    new Thread(() => {
        //        this.Dispatcher.Invoke(new Action(() => { statBarText.Text = "Exit the Application"; }));
        //    }).Start();
                      
        //}

        //protected void MouseLeaveArea(object sender, RoutedEventArgs args)
        //{
        //    statBarText.Text = "Ready";
        //}
        #endregion
       
        private void QuantificationRun_Click(object sender, RoutedEventArgs e)
        {
            if(!CheckParameters())
            {
                return;
            }

            System.DateTime data = System.DateTime.Now;
            string Now = data.ToString("yyyyMMdd-HH-mm-ss");
            string strParametersPath = txtResultPath.Text+"\\parameters"+Now+".params";
            if(!SaveParameters(strParametersPath))
            {
                return;
            }

            // So make sure that the identified proteins contain identifier; 
            try
            {
                string strIdentifiedProteinsPath;
                ComboBoxItem cbiIfContanStand = (ComboBoxContainStand.SelectedItem as ComboBoxItem);
                if (MaxQuantCheckBox.IsChecked == true&&cbiIfContanStand.Content.ToString()=="yes")
                {
                    strIdentifiedProteinsPath = txtInputDirectoryPath.Text + "\\proteinGroups.txt";
                    FileStream aFile = new FileStream(strIdentifiedProteinsPath, FileMode.Open);
                    StreamReader streamReader = new StreamReader(aFile);
                    ObservableCollection<ProteinName> ProteinNames = new ObservableCollection<ProteinName>();
                    string strLine;
                    int iBegain, iEnd;
                    string strProteinIdTemp;
                    bool bIdentifierIfExist = false;

                    strLine = streamReader.ReadLine();
                    strLine = streamReader.ReadLine();
                    while (strLine != null)
                    {
                        while (strLine == "")
                        {
                            strLine = streamReader.ReadLine();
                        }

                        iBegain = strLine.IndexOf("\t");
                        if (iBegain == -1)
                        {
                            System.Windows.MessageBox.Show("Cannot read protien Ids from proteinGroups.txt");
                        }
                        iEnd = strLine.IndexOf("\t", iBegain + 1);
                        if (iEnd == -1)
                        {
                            System.Windows.MessageBox.Show("Cannot read protien Ids from proteinGroups.txt");
                        }
                        strProteinIdTemp = strLine.Substring(iBegain + 1, iEnd - iBegain - 1);  // Assume that the second column is "Majority protein IDs"
                        if (strProteinIdTemp.IndexOf(txtIdentifierOfStand.Text) != -1)
                        {
                            bIdentifierIfExist = true;
                            break;
                        }

                        strLine = streamReader.ReadLine();
                    }

                    streamReader.Close();
                    if (!bIdentifierIfExist)
                    {
                        System.Windows.MessageBox.Show("Cannot find the standand protein identifier in input files.");
                        return;
                    }
                }
                else if(mzQuantMLCheckBox.IsChecked==true)
                {
                    // ToDo: Test if there is standard protein.
                }

            }
            catch (IOException ex)
            {
                System.Windows.MessageBox.Show("An IOException has been thrown!" + ex.Message + "\n Please set the input directory in advance!");
                return;
            }
            catch (ArgumentException ex)
            {
                System.Windows.MessageBox.Show("Cannot open " + ex.Message);
                return;
            }

            // set the parameter file for core program 
            string[] strParameters = new string[1];
            if((strParametersPath.IndexOf(" ")!=-1))
            {
                strParametersPath = "\"" + strParametersPath + "\"";
            }
            strParameters[0] = strParametersPath;            
            ///*************
            StartProcess(strParameters);     

        }
       
        private bool CheckParameters()
        {
            try
            {
                var inFiles = Directory.GetFiles(txtResultPath.Text);
            }
            catch (ArgumentException ex)
            {
                System.Windows.MessageBox.Show("The result directory cannot be empty!");
                return false;
            }
            catch (IOException ex)
            {
                System.Windows.MessageBox.Show("Cannot open the result directory: " + 
                    txtResultPath.Text);
                return false;
            }

            if(mzQuantMLCheckBox.IsChecked==false&&MaxQuantCheckBox.IsChecked==false
                &&SWATHCheckBox.IsChecked==false)
            {
                System.Windows.MessageBox.Show("Please choose one input type.");
                return true;
            }
            return true;
        }

        // return value: meaning
        // 0: success;
        // 1: Other errors;
        // 2: Cannot open the parameter file
        // 3: error when calling the C++ program;
        public void StartProcess( string[] args)
        {
            
            try
            {
                if (bIfQuantiRunning||bIfLoadRunning)
                {
                    return;
                }
               
                bIfExportResult = false;
                string s = "";
                foreach (string arg in args)
                {
                    s = s + arg + " ";
                }
                s = s.Trim();
                //new Thread(() =>
                //{
                //    //this.Dispatcher.BeginInvoke(new Action(() => { statBarText.Text = "Loading the input files"; }));
                //}).Start();

                //this.Dispatcher.BeginInvoke(System.Windows.Threading.DispatcherPriority.Normal,
                //new ChangeStatusBar(UpdateStatusBar), "Loading the input files");

                DisableParameters();
                LoadProcess=new Process();
                string loadExeName="Load.exe";

                ProcessStartInfo LoadStartInfo = new ProcessStartInfo(loadExeName, s);
                LoadStartInfo.CreateNoWindow = true;
                LoadStartInfo.UseShellExecute = false;
               //LoadStartInfo.RedirectStandardOutput = true;
                //hide  the window of exe
                //LoadStartInfo.WindowStyle = ProcessWindowStyle.Hidden;
                // run exe
                //Process proBach = Process.Start(LoadStartInfo);
                //proBach.WaitForExit();
               bIfRunLoad = true;

               UpdateStatusBar("Loading the input files");
               DoEvents();
               LoadProcess.EnableRaisingEvents = true;
              // LoadProcess.SynchronizingObject = this; 
               LoadProcess.Exited += new EventHandler(LoadProcess_Exited);
               LoadProcess.StartInfo = LoadStartInfo;
               bIfLoadRunning = true;
               LoadProcess.Start();
               LoadProcess.WaitForExit();
               
              
                if(LoadProcess.ExitCode!=0)
                {
                    EnableParameters();
                    return;
                }
                //hide  the window of exe
               

                string filename = "ProteinAbsoluteQuan.exe";
                ProcessStartInfo QuantiStartInfo = new ProcessStartInfo(filename, s);
                QuantiStartInfo.CreateNoWindow = true;
                QuantiStartInfo.UseShellExecute = false;
                //QuantiStartInfo.RedirectStandardOutput = true;
               // QuantiStartInfo.UseShellExecute = false;
                //QuantiStartInfo.WindowStyle = ProcessWindowStyle.Hidden;
                bIfRunQuant = true;
                QuantproBach = new Process();
                QuantproBach.EnableRaisingEvents = true;
                QuantproBach.Exited += null;
                QuantproBach.Exited += new EventHandler(QuantproBach_Exited);
                QuantproBach.StartInfo = QuantiStartInfo;
                bIfQuantiRunning = true;
              
                QuantproBach.Start();              
                statusBar.Background = System.Windows.Media.Brushes.Green;
                UpdateStatusBar("LFAQ is running! ");
                DoEvents();
                            
            }
            catch (Exception ex)
            {
                System.Windows.MessageBox.Show("Error when starting the program. The reason may be：" + ex.Message);
            }
        }

        private void LoadProcess_Exited(object sender, EventArgs e)
        {
            int returnLoadValue = LoadProcess.ExitCode;
            bIfLoadRunning = false;
            if (returnLoadValue == 2)
            {
                System.Windows.MessageBox.Show("Cannot open the parameter file! (The file path cannot contain space.)");               
            }
            else if (returnLoadValue != 0)
            {
                this.Dispatcher.BeginInvoke(System.Windows.Threading.DispatcherPriority.Normal,
new ChangeStatusBar(UpdateStatusBar), "Ready");
               // UpdateStatusBar("Ready");
                DoEvents();
                this.Dispatcher.BeginInvoke(System.Windows.Threading.DispatcherPriority.Normal,
new DelegateShowError(ShowError));                
            }
            else
            {
                //System.Windows.MessageBox.Show("Have loaded the input files!");
            }
            //throw new NotImplementedException();
        }

        private void QuantproBach_Exited(object sender, EventArgs e)
        {

           int returnLoadValue = QuantproBach.ExitCode;
           bIfQuantiRunning = false;
           this.Dispatcher.BeginInvoke(System.Windows.Threading.DispatcherPriority.Normal,
new ChangeParamIsEnable(EnableParameters));
           
            if (returnLoadValue != 0)
            {
                bIfQuanted = false;
                //UpdateStatusBar("Ready");
                this.Dispatcher.BeginInvoke(System.Windows.Threading.DispatcherPriority.Normal,
new ChangeStatusBar(UpdateStatusBar), "Ready");
                DoEvents();
                this.Dispatcher.BeginInvoke(System.Windows.Threading.DispatcherPriority.Normal,
new DelegateShowError(ShowError));
                
            }
            else
            {
                bIfQuanted = true;
                if(!bIfExportResult)
                {
                    
                    System.Windows.MessageBox.Show("LFAQ quantification finished successfully!");
                    bIfExportResult = true;
                    this.Dispatcher.BeginInvoke(System.Windows.Threading.DispatcherPriority.Normal,
new ChangeStatusBar(UpdateStatusBar), "Finished!");
                    //UpdateStatusBar("Finished! ");
                    DoEvents();
                }              

               
            }
            //throw new NotImplementedException();
        }
         public int KillProcesses()
        {
            try
            {
                if (bIfRunLoad)
                {
                    if(!LoadProcess.HasExited)
                        LoadProcess.Kill();
                  
                }
                if(bIfRunQuant)
                {
                    if(!QuantproBach.HasExited)
                        QuantproBach.Kill();
                }
            }
             catch(NullReferenceException nre)
            {

            }
                              
            return 0;
        }
        private void Parameters_Save_Click(object sender, RoutedEventArgs e)
        {
            System.Windows.Forms.SaveFileDialog ifile = new System.Windows.Forms.SaveFileDialog();
            ifile.Filter = "text files|*.params";
            if (ifile.ShowDialog() == System.Windows.Forms.DialogResult.OK)
            {
                SaveParameters(ifile.FileName);
            }
        }
        public bool SaveParameters( string path)
        {        
            using (StreamWriter writer = File.CreateText(path))
            {
                if(MaxQuantCheckBox.IsChecked==true)
                {
                    writer.WriteLine("IdentificationFileType=\"maxquant\"");
                    if (txtInputDirectoryPath.Text != "")
                    {
                        writer.WriteLine("Input directory path=\"" + txtInputDirectoryPath.Text + "\"");
                    }
                    else
                    {
                        System.Windows.MessageBox.Show("The input directory cannot be empty!");
                        return false;
                    }
                    try
                    {
                        string[] inFiles = Directory.GetFiles(txtInputDirectoryPath.Text);
                        List<string> ListFilesName = new List<string>();
                        int begain = 0;
                        foreach (string item in inFiles)
                        {
                            begain = item.LastIndexOf("\\");
                            ListFilesName.Add(item.Substring(begain + 1, item.Length - begain - 1));
                        }
                        if (!ListFilesName.Contains("peptides.txt"))
                        {
                            System.Windows.MessageBox.Show("Cannot find peptides.txt in the input directory.");
                            return false;
                        }
                        if (!ListFilesName.Contains("proteinGroups.txt"))
                        {
                            System.Windows.MessageBox.Show("Cannot find proteinGroups.txt in the input directory.");
                            return false;
                        }
                        if (!ListFilesName.Contains("experimentalDesignTemplate.txt"))
                        {
                            System.Windows.MessageBox.Show("Cannot find experimentalDesignTemplate.txt in the input directory.");
                            return false;
                        }
                    }
                    catch (ArgumentException ex)
                    {
                        System.Windows.MessageBox.Show("Cannot open the input directory: " + txtInputDirectoryPath.Text);
                        return false;
                    }
                    catch (IOException ex)
                    {
                        System.Windows.MessageBox.Show("Cannot open the input directory: " + txtInputDirectoryPath.Text);
                        return false;
                    }
                }
                else if (mzQuantMLCheckBox.IsChecked == true)
                {
                    writer.WriteLine("IdentificationFileType=\"mzQuantML\"");

                    if (!File.Exists(txtInputFilePath.Text))
                    {
                        System.Windows.MessageBox.Show("Cannot open the mzQuantML file: " + txtInputFilePath.Text);
                        return false;
                    }

                    if (txtInputFilePath.Text != "")
                    {
                        writer.WriteLine("Input file path=\"" + txtInputFilePath.Text + "\"");
                    }
                    else
                    {
                        System.Windows.MessageBox.Show("The path of mzQuantML file cannot be empty!");
                        return false;
                    }
                }
                else if (SWATHCheckBox.IsChecked==true)
                {
                    writer.WriteLine("IdentificationFileType=\"PeakView\"");

                    if (!File.Exists(txtInputFilePath.Text))
                    {
                        System.Windows.MessageBox.Show("Cannot open the SWATH file: " + txtInputFilePath.Text);
                        return false;
                    }

                    if (txtInputFilePath.Text != "")
                    {
                        writer.WriteLine("Input file path=\"" + txtInputFilePath.Text + "\"");
                    }
                    else
                    {
                        System.Windows.MessageBox.Show("The path of mzQuantML file cannot be empty!");
                        return false;
                    }

                }
                else
                {
                    System.Windows.MessageBox.Show("Please choose the indentification software at first!");
                    return false;
                }                    
                
                if (txtfasta.Text != "")
                {
                    writer.WriteLine("Fastapath=\"" + txtfasta.Text + "\"");
                }
                else
                {
                    System.Windows.MessageBox.Show("The identifier parsing rule cannot be empty!");
                    return false;
                }
                if (!File.Exists(txtfasta.Text))
                {
                    System.Windows.MessageBox.Show("Cannot open the fasta file: " + txtfasta.Text);
                    return false;
                }

                if (txtfastaType.Text != "")
                {
                    writer.WriteLine("IdentifierParsingRule=\"" + txtfastaType.Text + "\"");
                }
                else
                {
                    System.Windows.MessageBox.Show("The fasta Type cannot be empty!");
                    return false;
               }

                if(ContainDecoyProtein.IsChecked==true)
                {
                    writer.WriteLine("IfExistDecoyProteins=\"true\"");
                    if(txtDecoyPrefix.Text!="")
                    {
                        writer.WriteLine("PrefixOfDecoyProtein=\"" + txtDecoyPrefix.Text + "\"");
                    }
                    else
                    {
                        System.Windows.MessageBox.Show("The prefix of decoy proteins cannot be empty when the sample contains decoy proteins!");
                        return false;
                    }
                }
                else if(NoDecoyProtein.IsChecked==true)
                {
                    writer.WriteLine("IfExistDecoyProteins=\"false\"");
                    writer.WriteLine("PrefixOfDecoyProtein=\" \"");
                }
                else
                {
                    System.Windows.MessageBox.Show("Please choose if the fasta file contains decoy proteins or not!");
                    return false;
 
                }

                if (ExistContaProteins.IsChecked == true)
                {
                    writer.WriteLine("IfExistContaminantProteins=\"true\"");
                    if (txtContaminantPrefix.Text != "")
                    {
                        writer.WriteLine("PrefixOfContaminantProtein=\"" + txtContaminantPrefix.Text + "\"");
                    }
                    else
                    {
                        System.Windows.MessageBox.Show("The prefix of contaminant proteins cannot be empty when the sample contains contaminant proteins!");
                        return false;
                    }  
                }
                else if (NoContaProteins.IsChecked == true)
                {
                    writer.WriteLine("IfExistContaminantProteins=\"false\"");
                    writer.WriteLine("PrefixOfContaminantProtein=\" \"");
                }
                else
                {
                    System.Windows.MessageBox.Show("Please choose if the fasta file contains contaminant proteins or not!");
                    return false;
                }
                if(IfCalculateiBAQ.IsChecked==true)
                {
                    writer.WriteLine("IfCalculateiBAQ=\"true\"");
                }
                else
                {
                    writer.WriteLine("IfCalculateiBAQ=\"false\"");
                }
                if (IfCalculateTop3.IsChecked == true)
                {
                    writer.WriteLine("IfCalculateTop3=\"true\"");
                }
                else
                {
                    writer.WriteLine("IfCalculateTop3=\"false\"");
                }
                if (txtResultPath.Text != "")
                {
                    writer.WriteLine("ResultPath=\"" + txtResultPath.Text + "\"");
                }
                else
                {
                    System.Windows.MessageBox.Show("The result directory cannot be empty!");
                    return false;                  
                }
                
                ComboBoxItem cbi = (ComboBoxRegressionMethod.SelectedItem as ComboBoxItem);
                if (cbi.Content.ToString() == "BART")
                {
                    writer.WriteLine("RegressionMethod=\"BART\"");
                    if (txtAlpha.Text != "")
                    {
                        writer.WriteLine("alpha=\"" + txtAlpha.Text + "\"");
                    }
                    else
                    {
                        System.Windows.MessageBox.Show("The alpha cannot be empty!");
                        return false;
                    }
                    if (txtBeta.Text != "")
                    {
                        writer.WriteLine("beta=\"" + txtBeta.Text + "\"");
                    }
                    else
                    {
                        System.Windows.MessageBox.Show("The beta cannot be empty!");
                        return false;
                    }
                    if (txtK.Text != "")
                    {
                        writer.WriteLine("k=\"" + txtK.Text + "\"");
                    }
                    else
                    {
                        System.Windows.MessageBox.Show("The k cannot be empty!");
                        return false;
                    }
                    if (txtNumberOfTrees.Text != "")
                    {
                        writer.WriteLine("Number of trees=\"" + txtNumberOfTrees.Text + "\"");
                    }
                    else
                    {
                        System.Windows.MessageBox.Show("The Number of trees cannot be empty!");
                        return false;
                    }
                }
                else if (cbi.Content.ToString() == "stepwise")
                {
                    writer.WriteLine("RegressionMethod=\"stepwise\"");

                    if (txtAlpha1.Text != "")
                    {
                        writer.WriteLine("alpha1=\"" + txtAlpha1.Text + "\"");
                    }
                    else
                    {
                        System.Windows.MessageBox.Show("The alpha1 cannot be empty!");
                        return false;
                    }
                    if (txtAlpha2.Text != "")
                    {
                        writer.WriteLine("alpha2=\"" + txtAlpha2.Text + "\"");
                    }
                    else
                    {
                        System.Windows.MessageBox.Show("The alpha cannot be empty!");
                        return false;
                    }
                }

               
                if (txtMaxMissedCleave.Text != "")
                {
                    writer.WriteLine("MaxMissedCleavage=\"" + txtMaxMissedCleave.Text + "\"");
                }
                else
                {
                    System.Windows.MessageBox.Show("The maximun of missing cleavage cannot be empty!");
                    return false;
                }
                if (txtPepShotest.Text != "")
                {
                    writer.WriteLine("PepShortestLen=\"" + txtPepShotest.Text + "\"");
                }
                else
                {
                    System.Windows.MessageBox.Show("The shortest length of peptide cannot be empty!");
                    return false;
                }
                if (txtPepLongest.Text != "")
                {
                    writer.WriteLine("PepLongestLen=\"" + txtPepLongest.Text + "\"");
                }
                else
                {
                    System.Windows.MessageBox.Show("The longest length of peptide cannot be empty!");
                    return false;
                }
                if (txtEnzyme.Text != "")
                {
                    writer.WriteLine("Enzyme=\"" + txtEnzyme.Text + "\"");
                }
                else
                {
                    System.Windows.MessageBox.Show("The enzyme type cannot be empty!");
                    return false;
                }

                //if (CheckBoxOptimize.IsChecked==true)
                //{
                //    writer.WriteLine("IfOptimizeParameters=\"true\"");
                //}
                //else
                //{
                //    writer.WriteLine("IfOptimizeParameters=\"false\"");
                //}

                cbi = (ComboBoxContainStand.SelectedItem as ComboBoxItem);
                if (cbi.Content.ToString() == "yes")
                {
                    writer.WriteLine("IfCotainStandardProtein=\"true\"");
                    if (txtIdentifierOfStand.Text != "")
                    {
                        writer.WriteLine("IdentifierOfStandardProtein=\"" + txtIdentifierOfStand.Text + "\"");
                    }
                    else
                    {
                        System.Windows.MessageBox.Show("The identifiers of standard proteins cannot be empty!");
                        return false;
                    }
                }
                else if(cbi.Content.ToString()=="no")
                {
                    writer.WriteLine("IfCotainStandardProtein=\"false\"");

                }
                writer.WriteLine("StandardProteinsFilePath=\"" + txtResultPath.Text + "\\StandardProteins.txt\"");               
                writer.Close();
            }
            

            return true;
        }

        private void InputDirectory_Browse_Click(object sender, RoutedEventArgs e)
        {
          
            FolderBrowserDialog FBD = new FolderBrowserDialog();
            FBD.Description = "Please choose a directory:";
            if (FBD.ShowDialog() == System.Windows.Forms.DialogResult.OK)
            {
                txtInputDirectoryPath.Text = FBD.SelectedPath;
            }
        }
        private void fasta_browse_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openDlg = new Microsoft.Win32.OpenFileDialog();
            openDlg.Filter = "fasta files| *.fasta";
            if (true == openDlg.ShowDialog())
            {
                string DatafromFile = openDlg.FileName;
                txtfasta.Text = DatafromFile;
            }

        }
        
        private void fasta_type_Click(object sender, RoutedEventArgs e)
        {
            if(MaxQuantCheckBox.IsChecked==true)
            {
                try
                {
                    string[] inFiles = Directory.GetFiles(txtInputDirectoryPath.Text);
                    List<string> ListFilesName = new List<string>();
                    int begain = 0;
                    foreach (string item in inFiles)
                    {
                        begain = item.LastIndexOf("\\");
                        ListFilesName.Add(item.Substring(begain + 1, item.Length - begain - 1));
                    }
                    if (!ListFilesName.Contains("peptides.txt"))
                    {
                        System.Windows.MessageBox.Show("Cannot find peptides.txt in the input directory.");
                        return;
                    }
                    if (!ListFilesName.Contains("proteinGroups.txt"))
                    {
                        System.Windows.MessageBox.Show("Cannot find proteinGroups.txt in the input directory.");
                        return;
                    }
                    if (!ListFilesName.Contains("experimentalDesignTemplate.txt"))
                    {
                        System.Windows.MessageBox.Show("Cannot find experimentalDesignTemplate.txt in the input directory.");
                        return;
                    }
                }
                catch (ArgumentException ex)
                {
                    System.Windows.MessageBox.Show("Please set the input directory in advance!");
                    return;
                }
                catch (IOException ex)
                {
                    System.Windows.MessageBox.Show("Cannot open the input directory: " + txtInputDirectoryPath.Text + "Please set the input directory correctly!");
                    return;
                }
            }
            else if(SWATHCheckBox.IsChecked ==true)
            {
                if(txtInputFilePath.Text=="")
                {
                    System.Windows.MessageBox.Show("Please set the input file path in advance!");
                    return;
                }
                try
                {
                    if(!File.Exists(txtInputFilePath.Text))
                    {
                        System.Windows.MessageBox.Show("Please check the input file path again!");
                        return;
                    }                   
                }
                catch (ArgumentException ex)
                {
                    System.Windows.MessageBox.Show("Please set the input file path in advance!");
                    return;
                }
                catch (IOException ex)
                {
                    System.Windows.MessageBox.Show("Cannot open the input file: " + txtInputDirectoryPath.Text + "Please set the input file path correctly!");
                    return;
                }
            }
            
            

            SetRegularExpression re = new SetRegularExpression();
            re.ChangeTextEvent += new ChangeTextHandler(SetCorrectRE);
            re.getInputPathText += new GetTextBoxContentHandler(GettxtInputPath);
            re.getInputType += new GetInputTypeHandle(GetInputType);
            re.Owner = this;
            re.ShowDialog();
        }
        private void SetStandProteins_Click(object sender, RoutedEventArgs e)
        {
            if(MaxQuantCheckBox.IsChecked==true)
            {
                try
                {
                    string[] inFiles = Directory.GetFiles(txtInputDirectoryPath.Text);
                    List<string> ListFilesName = new List<string>();
                    int begain = 0;
                    foreach (string item in inFiles)
                    {
                        begain = item.LastIndexOf("\\");
                        ListFilesName.Add(item.Substring(begain + 1, item.Length - begain - 1));
                    }
                    if (!ListFilesName.Contains("peptides.txt"))
                    {
                        System.Windows.MessageBox.Show("Cannot find peptides.txt in the input directory.");
                        return;
                    }
                    if (!ListFilesName.Contains("proteinGroups.txt"))
                    {
                        System.Windows.MessageBox.Show("Cannot find proteinGroups.txt in the input directory.");
                        return;
                    }
                    if (!ListFilesName.Contains("experimentalDesignTemplate.txt"))
                    {
                        System.Windows.MessageBox.Show("Cannot find experimentalDesignTemplate.txt in the input directory.");
                        return;
                    }
                }
                catch (ArgumentException ex)
                {
                    System.Windows.MessageBox.Show("Please set the input directory in advance!");
                    return;
                }
                catch (IOException ex)
                {
                    System.Windows.MessageBox.Show("Cannot open the input directory: " + txtInputDirectoryPath.Text + "Please set the input directory correctly!");
                    return;
                }
            }
            
             
            if(txtIdentifierOfStand.Text=="")
            {
                System.Windows.MessageBox.Show("Please set the identifier of standard proteins in advance!");
                return;
            }

            SetStandProteins WinsetStand = new SetStandProteins();
            WinsetStand.GetIdentifierOfStand += new GetTextBoxContentHandler(GetIdentifierOfStand);
            WinsetStand.GetInputPath += new GetInputPathHandler(GettxtInputPath);
            WinsetStand.GetInputType += new GetInputTypeHandle(GetInputType);
            WinsetStand.GetResultPath += new GetResultPathHandler(GetResultPath);
            WinsetStand.Owner = this;
            
            WinsetStand.ShowDialog();
            //WinsetStand.Show();
        }
        void SetCorrectRE(string text)
        {
            this.txtfastaType.Text = text;
        }
        private void Result_Browse_Click(object sender, RoutedEventArgs e)
        {
            FolderBrowserDialog FBD = new FolderBrowserDialog();
            FBD.Description = "Please choose a directory: ";
            if (FBD.ShowDialog() == System.Windows.Forms.DialogResult.OK)
            {
                txtResultPath.Text = FBD.SelectedPath;
            }
        }


        private void method_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            ComboBoxItem cbi = ((sender as System.Windows.Controls.ComboBox).SelectedItem as ComboBoxItem);

            if(cbi.Content.ToString()=="BART")
            {
                if (stepwisePanel != null)
                {
                    stepwisePanel.Visibility = Visibility.Hidden;
                }
                if(BARTPanel!=null)
                {
                    BARTPanel.Visibility=Visibility.Visible;
                }
            }
            else if(cbi.Content.ToString()=="stepwise")
            {
                if (BARTPanel != null)
                {
                    BARTPanel.Visibility = Visibility.Hidden;
                }
                if (stepwisePanel != null)
                {
                    stepwisePanel.Visibility = Visibility.Visible;
                }

            }  
        }

        public void LoadParameters(string parametersPath)
        {
            string strLine;
            string strTemp1, strTemp2;
            int iSplitLoc;
            int iStartLoc, iEndLoc;
            try               
            {
                StreamReader sr = File.OpenText(parametersPath);
                while (sr.EndOfStream != true)
                {
                    strLine = sr.ReadLine();
                    iSplitLoc = strLine.IndexOf("=");
                    if (iSplitLoc == -1)
                    {
                        continue;
                    }
                    strTemp1 = strLine.Substring(0, iSplitLoc);
                    iStartLoc = strLine.IndexOf("\"");
                    iEndLoc = strLine.LastIndexOf("\"");
                    strTemp2 = strLine.Substring(iStartLoc + 1, iEndLoc - iStartLoc - 1);
                    if (strTemp1 == "IdentificationFileType")
                    {
                        if (strTemp2 == "maxquant")
                        {
                            MaxQuantCheckBox.IsChecked = true;
                        }
                        else if (strTemp2 == "mzQuantML")
                        {
                            mzQuantMLCheckBox.IsChecked = true;
                        }
                        else if(strTemp2=="PeakView")
                        {
                            SWATHCheckBox.IsChecked = true;
                        }

                    }
                    if (strTemp1 == "Input directory path")
                    {
                        txtInputDirectoryPath.Text = strTemp2;
                    }
                    if (strTemp1 == "Input file path")
                    {
                        txtInputFilePath.Text = strTemp2;
                    }
                    if (strTemp1 == "Fastapath")
                    {
                        txtfasta.Text = strTemp2;
                    }
                    if (strTemp1 == "IdentifierParsingRule")
                    {
                        txtfastaType.Text = strTemp2;
                    }

                    if (strTemp1 == "IfExistDecoyProteins")
                    {
                        if (strTemp2 == "true")
                        {
                            ContainDecoyProtein.IsChecked = true;
                        }
                        else if (strTemp2 == "false")
                        {
                            NoDecoyProtein.IsChecked = true;
                        }
                    }
                    if (strTemp1 == "PrefixOfDecoyProtein")
                    {
                        txtDecoyPrefix.Text = strTemp2;
                    }
                    if (strTemp1 == "IfExistContaminantProteins")
                    {
                        if (strTemp2 == "true")
                        {
                            ExistContaProteins.IsChecked = true;
                        }
                        else if (strTemp2 == "false")
                        {
                            NoContaProteins.IsChecked = true;
                        }
                    }
                    if (strTemp1 == "PrefixOfContaminantProtein")
                    {
                        txtContaminantPrefix.Text = strTemp2;
                    }
                    if (strTemp1 == "IfCalculateiBAQ")
                    {
                        if (strTemp2 == "true")
                        {
                            IfCalculateiBAQ.IsChecked = true;
                        }
                        else if (strTemp2 == "false")
                        {
                            IfCalculateiBAQ.IsChecked = true;
                        }
                    }
                    if (strTemp1 == "IfCalculateTop3")
                    {
                        if (strTemp2 == "true")
                        {
                            IfCalculateTop3.IsChecked = true;
                        }
                        else if (strTemp2 == "false")
                        {
                            IfCalculateTop3.IsChecked = true;
                        }
                    }

                    if (strTemp1 == "RegressionMethod")
                    {
                        if (strTemp2 == "BART")
                        {
                            ComboBoxRegressionMethod.SelectedIndex = 1;
                        }
                        else if (strTemp2 == "stepwise")
                        {
                            ComboBoxRegressionMethod.SelectedIndex = 0;
                        }
                    }

                    if (strTemp1 == "alpha1")
                    {
                        txtAlpha1.Text = strTemp2;
                    }
                    if (strTemp1 == "Alpha2")
                    {
                        txtAlpha2.Text = strTemp2;
                    }
                    if (strTemp1 == "alpha")
                    {
                        txtAlpha.Text = strTemp2;
                    }

                    if (strTemp1 == "beta")
                    {
                        txtBeta.Text = strTemp2;
                    }
                    if (strTemp1 == "Number of trees")
                    {
                        txtNumberOfTrees.Text = strTemp2;
                    }
                    if (strTemp1 == "k")
                    {
                        txtK.Text = strTemp2;
                    }
                    if (strTemp1 == "MaxMissedCleavage")
                    {
                        txtMaxMissedCleave.Text = strTemp2;
                    }

                    if (strTemp1 == "PepShortestLen")
                    {
                        txtPepShotest.Text = strTemp2;
                    }
                    if (strTemp1 == "PepLongestLen")
                    {
                        txtPepLongest.Text = strTemp2;
                    }

                    if (strTemp1 == "ResultPath")
                    {
                        txtResultPath.Text = strTemp2;
                    }
                    //if (strTemp1 == "IfOptimizeParameters" && strTemp2 == "true")
                    //{

                    //    //CheckBoxOptimize.IsChecked = true;
                    //}
                    if (strTemp1 == "IfCotainStandardProtein")
                    {
                        if (strTemp2 == "true")
                        {
                            ComboBoxContainStand.SelectedIndex = 0;
                        }
                        else if (strTemp2 == "false")
                        {
                            ComboBoxContainStand.SelectedIndex = 1;
                        }

                    }
                    if (strTemp1 == "IdentifierOfStandardProtein")
                    {
                        txtIdentifierOfStand.Text = strTemp2;
                    }
                }
                sr.Close();
            }
            catch(IOException ex)
            {
                System.Windows.MessageBox.Show("Cannot open " + parametersPath);
            }
            catch (ArgumentException ax)
            {
                System.Windows.MessageBox.Show("Cannot open " + parametersPath);
            }

        }
       private void Load_paramters_Click(object sender, RoutedEventArgs e)
       {
          System.Windows.Forms.OpenFileDialog ofile = new System.Windows.Forms.OpenFileDialog();
           ofile.Filter = "textfiles|*.params";
           if(ofile.ShowDialog()==System.Windows.Forms.DialogResult.OK)
           {
               LoadParameters(ofile.FileName);
           }          
       }

       private void DecoyProteinChecked(object sender, RoutedEventArgs e)
       {
           try
           {
               NoDecoyProtein.IsChecked = false;
               txtDecoyPrefix.IsEnabled = true;
           }
           catch (IOException ex)
           {
             System.Windows.MessageBox.Show("An IOException has been thrown!"+ex.Message);
              return;
           }
       }

       private void DecoyProteinUnchecked(object sender, RoutedEventArgs e)
       {
           try
           {
               NoDecoyProtein.IsChecked = true;
               txtDecoyPrefix.IsEnabled = false;
           }
           catch (IOException ex)
           {
               System.Windows.MessageBox.Show("An IOException has been thrown!"+ex.Message);
               return;
           }

       }
       private void NoDecoyProteinChecked(object sender, RoutedEventArgs e)
       {
           try
           {
               ContainDecoyProtein.IsChecked = false;
               if(this.IsLoaded)
                    txtDecoyPrefix.IsEnabled = false;
           }
           catch (IOException ex)
           {
               System.Windows.MessageBox.Show("An IOException has been thrown!"+ex.Message);
               return;
           }
       }
       private void NoDecoyProteinUnchecked(object sender, RoutedEventArgs e)
       {
           try
           {
               ContainDecoyProtein.IsChecked = true;
               txtDecoyPrefix.IsEnabled = true;
           }
           catch (IOException ex)
           {
               System.Windows.MessageBox.Show("An IOException has been thrown!"+ex.Message);
               return;
           }

       }
       private void NoContaProteins_Checked(object sender, RoutedEventArgs e)
       {
           try
           {
               ExistContaProteins.IsChecked = false;
               if (this.IsLoaded)
                    txtContaminantPrefix.IsEnabled = false;
           }
           catch (IOException ex)
           {
               System.Windows.MessageBox.Show("An IOException has been thrown!"+ex.Message);
               return;
           }
       }
       private void NoContaProteins_Unchecked(object sender, RoutedEventArgs e)
       {
           try
           {
               ExistContaProteins.IsChecked = true;
               txtContaminantPrefix.IsEnabled = true;
           }
           catch (IOException ex)
           {
               System.Windows.MessageBox.Show("An IOException has been thrown!"+ex.Message);
               return;
           }
       }

       private void ExistContaProteinsUnchecked(object sender, RoutedEventArgs e)
       {
           try
           {
               NoContaProteins.IsChecked = true;
               txtContaminantPrefix.IsEnabled = false;
           }
           catch (IOException ex)
           {
               System.Windows.MessageBox.Show("An IOException has been thrown!"+ex.Message);
               return;
           }
       }
       private void ExistContaProteinsChecked(object sender, RoutedEventArgs e)
       {
           try
           {
               NoContaProteins.IsChecked = false;
               txtContaminantPrefix.IsEnabled = true;
           }
           catch (IOException ex)
           {
               System.Windows.MessageBox.Show("An IOException has been thrown!"+ex.Message);
               return;
           }
           
       }


        private void txtfasta_MouseLeave(object sender, System.Windows.Input.MouseEventArgs e)
        {
            string str = txtfasta.Text;
            if (str.Length > 200)
            {
                System.Windows.MessageBox.Show("The length of directory cannot be more than 200.");
            }
        }

        private void FilePath_MouseLeave(object sender, System.Windows.Input.MouseEventArgs e)
        {
            string str = txtInputFilePath.Text;
            if (str.Length > 200)
            {
                System.Windows.MessageBox.Show("The length of directory cannot be more than 200.");
            }
        }

        private void txtDirectoryPath_MouseLeave(object sender, System.Windows.Input.MouseEventArgs e)
        {
            string str = txtInputDirectoryPath.Text;
            if (str.Length > 200)
            {
                System.Windows.MessageBox.Show("The length of directory cannot be more than 200.");
            }
        }
        private void txtResultPath_MouseLeave(object sender, System.Windows.Input.MouseEventArgs e)
        {
            string str = txtResultPath.Text;
            if (str.Length > 200)
            {
                System.Windows.MessageBox.Show("The length of result directory cannot be more than 200.");
            }
        }

        private void txtAlpha_TextChanged(object sender, TextChangedEventArgs e)
        {
            string strTemp = txtAlpha.Text;
            if(!string.IsNullOrEmpty(strTemp))
            {
                try
                {
                    double dTemp = Convert.ToDouble(strTemp);
                    if (dTemp < 0.0 || dTemp > 1.0)
                    {
                        System.Windows.MessageBox.Show("The alpha should be in the 0-1 range");
                        e.Handled = false;
                    }
                }
                catch
                {
                    System.Windows.MessageBox.Show("The alpha should be a double.");
                }
            }
           
        }
        private void txtBeta_TextChanged(object sender, TextChangedEventArgs e)
        {

            string strTemp = txtAlpha.Text;
            if (!string.IsNullOrEmpty(strTemp))
            {
                try
                {
                    double dTemp = Convert.ToDouble(strTemp);
                    if (dTemp < 0.0 || dTemp > 3.0)
                    {
                        System.Windows.MessageBox.Show("The beta should be in the 0.0-3.0 range");
                        e.Handled = false;
                    }
                }
                catch
                {
                    System.Windows.MessageBox.Show("The beta should be a double.");
                }
            }
        }
        private void txtK_TextChanged(object sender, TextChangedEventArgs e)
        {

            string strTemp = txtK.Text;
            if (!string.IsNullOrEmpty(strTemp))
            {
                try
                {
                    Int32 dTemp = Convert.ToInt32(strTemp);
                    if (dTemp < 1 || dTemp > 5)
                    {
                        System.Windows.MessageBox.Show("The k should be in the 1-5 range");
                        e.Handled = false;
                    }
                }
                catch
                {
                    System.Windows.MessageBox.Show("The k should be a integer.");
                }
            }
        }
        private void txtNumberOfTrees_TextChanged(object sender, TextChangedEventArgs e)
        {

            string strTemp = txtNumberOfTrees.Text;
            if (!string.IsNullOrEmpty(strTemp))
            {
                try
                {
                    Int32 dTemp = Convert.ToInt32(strTemp);
                    if (dTemp < 1 || dTemp > 500)
                    {
                        System.Windows.MessageBox.Show("The number of trees should be in the 1-500 range");
                        e.Handled = false;
                    }
                }
                catch
                {
                    System.Windows.MessageBox.Show("The number of trees should be a integer.");
                }
            }
        }
        private void txtMaxMissedCleave_TextChanged(object sender, TextChangedEventArgs e)
        {

            string strTemp = txtMaxMissedCleave.Text;
            if (!string.IsNullOrEmpty(strTemp))
            {
                try
                {
                    Int32 dTemp = Convert.ToInt32(strTemp);
                    if (dTemp < 0 || dTemp > 5)
                    {
                        System.Windows.MessageBox.Show("The allowed maximum missing cleavages of one peptide in theoretic digestion should be in the 0-5 range");
                        e.Handled = false;
                    }
                }
                catch
                {
                    System.Windows.MessageBox.Show("The allowed maximum missing cleavages of one peptide in theoretic digestion should be a integer.");
                }
            }
        }
        private void txtPepShotest_TextChanged(object sender, TextChangedEventArgs e)
        {

            string strTemp = txtPepShotest.Text;
            if (!string.IsNullOrEmpty(strTemp))
            {
                try
                {
                    Int32 dTemp = Convert.ToInt32(strTemp);
                    if (dTemp < 1 || dTemp > 10)
                    {
                        System.Windows.MessageBox.Show("The allowed shortest length of one peptide in theoretic digestion should be in the 1-10 range");
                        e.Handled = false;
                    }
                }
                catch
                {
                    System.Windows.MessageBox.Show("The allowed shortest length of one peptide in theoretic digestion should be a integer.");
                }
            }
        }
        private void txtPepLongest_TextChanged(object sender, TextChangedEventArgs e)
        {

            string strTemp = txtPepLongest.Text;
            if (!string.IsNullOrEmpty(strTemp))
            {
                try
                {
                    Int32 dTemp = Convert.ToInt32(strTemp);
                    if (dTemp < 10 || dTemp > 100)
                    {
                        System.Windows.MessageBox.Show("The allowed longest length of one peptide in theoretic digestion should be in the 10-100 range");
                        e.Handled = false;
                    }
                }
                catch
                {
                    System.Windows.MessageBox.Show("The allowed longest length of one peptide in theoretic digestion should be a integer.");
                }
            }
        }
        private void txtAlpha1_TextChanged(object sender, TextChangedEventArgs e)
        {
            string strTemp = txtAlpha1.Text;
            if (!string.IsNullOrEmpty(strTemp))
            {
                try
                {
                    double dTemp = Convert.ToDouble(strTemp);
                    if (dTemp < 0.0 || dTemp > 1.0)
                    {
                        System.Windows.MessageBox.Show("The alpha1 should be in the 0-1 range");
                        e.Handled = false;
                    }
                }
                catch
                {
                    System.Windows.MessageBox.Show("The alpha1 should be a double.");
                }
            }
        }
        private void txtAlpha2_TextChanged(object sender, TextChangedEventArgs e)
        {
            string strTemp = txtAlpha2.Text;
            if (!string.IsNullOrEmpty(strTemp))
            {
                try
                {
                    double dTemp = Convert.ToDouble(strTemp);
                    if (dTemp < 0.0 || dTemp > 1.0)
                    {
                        System.Windows.MessageBox.Show("The alpha2 should be in the 0-1 range");
                        e.Handled = false;
                    }
                }
                catch
                {
                    System.Windows.MessageBox.Show("The alpha2 should be a double.");
                }
            }
        }

        private void txtAlpha_PreviewKeyDown(object sender, System.Windows.Input.KeyEventArgs e)
        {
            
            int PPos = this.txtAlpha.Text.IndexOf('.');//the location of the decimal
            if ((e.Key >= Key.D1 && e.Key <= Key.D9)//keyboard 1-9
                || (e.Key >= Key.NumPad1 && e.Key <= Key.NumPad9))//keypad 1-9
            {
                if (this.txtAlpha.SelectionStart > 0)
                {
                    if (this.txtAlpha.Text.Substring(0, 1) == "0" && PPos != 1)
                    {//if the first charater is 0 ,the second chracter can only be decimal
                        this.txtAlpha.Text = this.txtAlpha.Text.Remove(0, 1);
                    }
                }
                e.Handled = false;
                return;
            }
            else if (e.Key == Key.D0 || e.Key == Key.NumPad0)//0
            {
                if (this.txtAlpha.SelectionStart > 0)
                {//Not first key
                    if (this.txtAlpha.Text.Substring(0, 1) != "0")
                    {//First digit is not 0
                        e.Handled = false;
                        return;
                    }
                    else
                    {//First digit is 0
                        if (PPos == 1)
                        { //0.
                            e.Handled = false;
                            return;
                        }
                    }
                }
                else
                {
                    e.Handled = false;
                    return;
                }
            }
            else if (e.Key == Key.Decimal || e.Key == Key.OemPeriod)//decimal
            {
                if (PPos == -1)//there has not yet been a decimal
                {
                    if (this.txtAlpha.SelectionStart == 0)
                    {
                        this.txtAlpha.Text = "0.";
                        this.txtAlpha.Select(this.txtAlpha.Text.Length, 0);
                    }
                    e.Handled = false;
                    return;
                }
            }
            else if (e.Key == Key.Back)
            {
                e.Handled = false;
                return;
            }
            e.Handled = true;
        }
        private void txtBeta_PreviewKeyDown(object sender, System.Windows.Input.KeyEventArgs e)
        {
            int PPos = this.txtBeta.Text.IndexOf('.');//the location of the decimal
            if ((e.Key >= Key.D1 && e.Key <= Key.D9)//keyboard 1-9
                || (e.Key >= Key.NumPad1 && e.Key <= Key.NumPad9))//keypad 1-9
            {
                if (this.txtBeta.SelectionStart > 0)
                {
                    if (this.txtBeta.Text.Substring(0, 1) == "0" && PPos != 1)
                    {//if the first charater is 0 ,the second chracter can only be decimal
                        this.txtBeta.Text = this.txtBeta.Text.Remove(0, 1);
                    }
                }
                e.Handled = false;
                return;
            }
            else if (e.Key == Key.D0 || e.Key == Key.NumPad0)//0
            {
                if (this.txtBeta.SelectionStart > 0)
                {//Not first key
                    if (this.txtBeta.Text.Substring(0, 1) != "0")
                    {//First digit is not 0
                        e.Handled = false;
                        return;
                    }
                    else
                    {//First digit is 0
                        if (PPos == 1)
                        { //0.
                            e.Handled = false;
                            return;
                        }
                    }
                }
                else
                {
                    e.Handled = false;
                    return;
                }
            }
            else if (e.Key == Key.Decimal || e.Key == Key.OemPeriod)//decimal
            {
                if (PPos == -1)//there has not yet been a decimal
                {
                    if (this.txtBeta.SelectionStart == 0)
                    {
                        this.txtBeta.Text = "0.";
                        this.txtBeta.Select(this.txtBeta.Text.Length, 0);
                    }
                    e.Handled = false;
                    return;
                }
            }
            else if (e.Key == Key.Back)
            {
                e.Handled = false;
                return;
            }
            e.Handled = true;
        }
        private void txtK_PreviewKeyDown(object sender, System.Windows.Input.KeyEventArgs e)
        {
            if ((e.Key >= Key.D1 && e.Key <= Key.D9)//keyboard 1-9
                  || (e.Key >= Key.NumPad1 && e.Key <= Key.NumPad9))//keypad 1-9
            {
                e.Handled = false;
                return;
            }
            else if (e.Key == Key.Back)
            {
                e.Handled = false;
                return;
            }
            e.Handled = true;
        }
        private void txtNumberOfTrees_PreviewKeyDown(object sender, System.Windows.Input.KeyEventArgs e)
        {
            if ((e.Key >= Key.D1 && e.Key <= Key.D9)//keyboard 1-9
                 || (e.Key >= Key.NumPad1 && e.Key <= Key.NumPad9))//keypad 1-9
            {
                e.Handled = false;
                return;
            }
            else if (e.Key == Key.D0 || e.Key == Key.NumPad0)//0
            {
                if (this.txtNumberOfTrees.SelectionStart > 0)
                {//Not first key
                    if (this.txtNumberOfTrees.Text.Substring(0, 1) != "0")
                    {//First digit is not 0
                        e.Handled = false;
                        return;
                    }
                }
                else
                {
                    e.Handled = false;
                    return;
                }
            }
            else if (e.Key == Key.Back)
            {
                e.Handled = false;
                return;
            }
            e.Handled = true;
        }
        private void txtMaxMissedCleave_PreviewKeyDown(object sender, System.Windows.Input.KeyEventArgs e)
        {
            if ((e.Key >= Key.D0 && e.Key <= Key.D9)//keyboard 1-9
                 || (e.Key >= Key.NumPad1 && e.Key <= Key.NumPad9))//keypad 1-9
            {
                e.Handled = false;
                return;
            }
            else if (e.Key == Key.Back)
            {
                e.Handled = false;
                return;
            }
            e.Handled = true;
        }
        private void txtPepShotest_PreviewKeyDown(object sender, System.Windows.Input.KeyEventArgs e)
        {
            if ((e.Key >= Key.D1 && e.Key <= Key.D9)//keyboard 1-9
                 || (e.Key >= Key.NumPad1 && e.Key <= Key.NumPad9))//keypad 1-9
            {
                e.Handled = false;
                return;
            }
            else if (e.Key == Key.D0 || e.Key == Key.NumPad0)//0
            {
                if (this.txtPepShotest.SelectionStart > 0)
                {//Not first key
                    if (this.txtPepShotest.Text.Substring(0, 1) != "0")
                    {//First digit is not 0
                        e.Handled = false;
                        return;
                    }
                }
                else
                {
                    e.Handled = false;
                    return;
                }
            }
            else if (e.Key == Key.Back)
            {
                e.Handled = false;
                return;
            }
            e.Handled = true;
        }
        private void txtPepLongest_PreviewKeyDown(object sender, System.Windows.Input.KeyEventArgs e)
        {
            if ((e.Key >= Key.D1 && e.Key <= Key.D9)//keyboard 1-9
                 || (e.Key >= Key.NumPad1 && e.Key <= Key.NumPad9))//keypad 1-9
            {
                e.Handled = false;
                return;
            }
            else if (e.Key == Key.D0 || e.Key == Key.NumPad0)//0
            {
                if (this.txtPepLongest.SelectionStart > 0)
                {//Not first key
                    if (this.txtPepLongest.Text.Substring(0, 1) != "0")
                    {//First digit is not 0
                        e.Handled = false;
                        return;
                    }
                }
                else
                {
                    e.Handled = false;
                    return;
                }
            }
            else if (e.Key == Key.Back)
            {
                e.Handled = false;
                return;
            }
            e.Handled = true;
        }
        private void txtAlpha1_PreviewKeyDown(object sender, System.Windows.Input.KeyEventArgs e)
        {
            int PPos = this.txtAlpha1.Text.IndexOf('.');//the location of the decimal
            if ((e.Key >= Key.D1 && e.Key <= Key.D9)//keyboard 1-9
                || (e.Key >= Key.NumPad1 && e.Key <= Key.NumPad9))//keypad 1-9
            {
                if (this.txtAlpha1.SelectionStart > 0)
                {
                    if (this.txtAlpha1.Text.Substring(0, 1) == "0" && PPos != 1)
                    {//if the first charater is 0 ,the second chracter can only be decimal
                        this.txtAlpha1.Text = this.txtAlpha1.Text.Remove(0, 1);
                    }
                }
                e.Handled = false;
                return;
            }
            else if (e.Key == Key.D0 || e.Key == Key.NumPad0)//0
            {
                if (this.txtAlpha1.SelectionStart > 0)
                {//Not first key
                    if (this.txtAlpha1.Text.Substring(0, 1) != "0")
                    {//First digit is not 0
                        e.Handled = false;
                        return;
                    }
                    else
                    {//First digit is 0
                        if (PPos == 1)
                        { //0.
                            e.Handled = false;
                            return;
                        }
                    }
                }
                else
                {
                    e.Handled = false;
                    return;
                }
            }
            else if (e.Key == Key.Decimal || e.Key == Key.OemPeriod)//decimal
            {
                if (PPos == -1)//there has not yet been a decimal
                {
                    if (this.txtAlpha1.SelectionStart == 0)
                    {
                        this.txtAlpha1.Text = "0.";
                        this.txtAlpha1.Select(this.txtAlpha1.Text.Length, 0);
                    }
                    e.Handled = false;
                    return;
                }
            }
            else if (e.Key == Key.Back)
            {
                e.Handled = false;
                return;
            }
            e.Handled = true;
        }
        private void txtAlpha2_PreviewKeyDown(object sender, System.Windows.Input.KeyEventArgs e)
        {
            int PPos = this.txtAlpha2.Text.IndexOf('.');//the location of the decimal
            if ((e.Key >= Key.D1 && e.Key <= Key.D9)//keyboard 1-9
                || (e.Key >= Key.NumPad1 && e.Key <= Key.NumPad9))//keypad 1-9
            {
                if (this.txtAlpha2.SelectionStart > 0)
                {
                    if (this.txtAlpha2.Text.Substring(0, 1) == "0" && PPos != 1)
                    {//if the first charater is 0 ,the second chracter can only be decimal
                        this.txtAlpha2.Text = this.txtAlpha2.Text.Remove(0, 1);
                    }
                }
                e.Handled = false;
                return;
            }
            else if (e.Key == Key.D0 || e.Key == Key.NumPad0)//0
            {
                if (this.txtAlpha2.SelectionStart > 0)
                {//Not first key
                    if (this.txtAlpha2.Text.Substring(0, 1) != "0")
                    {//First digit is not 0
                        e.Handled = false;
                        return;
                    }
                    else
                    {//First digit is 0
                        if (PPos == 1)
                        { //0.
                            e.Handled = false;
                            return;
                        }
                    }
                }
                else
                {
                    e.Handled = false;
                    return;
                }
            }
            else if (e.Key == Key.Decimal || e.Key == Key.OemPeriod)//decimal
            {
                if (PPos == -1)//there has not yet been a decimal
                {
                    if (this.txtAlpha2.SelectionStart == 0)
                    {
                        this.txtAlpha2.Text = "0.";
                        this.txtAlpha2.Select(this.txtAlpha2.Text.Length, 0);
                    }
                    e.Handled = false;
                    return;
                }
            }
            else if (e.Key == Key.Back)
            {
                e.Handled = false;
                return;
            }
            e.Handled = true;
        }
        private void ParseHeadLine(string []strValues, out int iProteinIdColumnNum, out int iMajorProteinIdColumnNum,
            out int iExperimentColumnNum, out int iiBAQColumnNum, out int iLFAQpepColumnNum, 
            out int iTop3ColumnNum,out int iPredictMolOfLFAQColumnNum,out int iPredictMolOfiBAQColumnNum,
            out int iPredictMolOfTop3ColumnNum,out int iPeptidesIntensitiesColumnNum)
        {
           bool bIfFindProteinIdColumn=false; 
           bool bIfFindMajorProteinIdColumn=false;
           bool bIfFindExperimentColumn=false; 
           bool bIfFindiBAQColumn=false; 
           bool bIfFindLFAQpepColumn=false; 
           bool bIfFindTop3Column=false;
           bool bIfFindPredictMolOfLFAQColumn=false;
           bool bIfFindPredictMolOfiBAQColumn = false;
           bool bIfFindPredictMolOfTop3Column = false;
           bool bIfFindPeptidesIntensitiesColumn = false;

           iProteinIdColumnNum=0;
           iMajorProteinIdColumnNum=0;
           iExperimentColumnNum=0;
           iiBAQColumnNum=0;
           iLFAQpepColumnNum=0; 
           iTop3ColumnNum=0;
           iPredictMolOfLFAQColumnNum=0;
           iPredictMolOfiBAQColumnNum=0;
           iPredictMolOfTop3ColumnNum=0;
           iPeptidesIntensitiesColumnNum = 0;

           int indexTemp=0;
            foreach(string str in strValues)
            {
                if(str=="Protein IDs")
                {
                    iProteinIdColumnNum = indexTemp;
                    bIfFindProteinIdColumn = true;
                }
                else if(str=="Majority protein IDs")
                {
                    iMajorProteinIdColumnNum = indexTemp;
                    bIfFindMajorProteinIdColumn = true;
                }
                else if(str=="Experiments")
                {
                    iExperimentColumnNum = indexTemp;
                    bIfFindExperimentColumn = true;
                }
                else if(str=="LFAQ")
                {
                    iLFAQpepColumnNum = indexTemp;
                    bIfFindLFAQpepColumn = true;
                }
                else if(str=="iBAQ")
                {
                    iiBAQColumnNum = indexTemp;
                    bIfFindiBAQColumn = true;
                }
                else if(str=="Top3")
                {
                    iTop3ColumnNum = indexTemp;
                    bIfFindTop3Column = true;
                }
                else if(str=="PredictedMol(LFAQ)")
                {
                    iPredictMolOfLFAQColumnNum = indexTemp;
                    bIfFindPredictMolOfLFAQColumn = true;
                }
                else if(str=="PredictedMol(iBAQ)")
                {
                    iPredictMolOfiBAQColumnNum = indexTemp;
                    bIfFindPredictMolOfiBAQColumn = true;
                }
                else if(str=="PredictedMol(Top3)")
                {
                    iPredictMolOfTop3ColumnNum = indexTemp;
                    bIfFindPredictMolOfTop3Column = true;
                }
                else if(str=="PeptidesIntensities")
                {
                    iPeptidesIntensitiesColumnNum = indexTemp;
                    bIfFindPeptidesIntensitiesColumn = true;
                }
                indexTemp++;
            }

           if(bIfFindProteinIdColumn == false)
           {
               System.Windows.MessageBox.Show("Cannot find column \"Protein IDs\" in the ProteinMergedResults.txt\n");
               return;
           }

          if(bIfFindMajorProteinIdColumn == false)
          {
              System.Windows.MessageBox.Show("Cannot find column \"Majority protein IDs\" in the ProteinMergedResults.txt\n");
              return;
          }
          if(bIfFindExperimentColumn == false)
          {
              System.Windows.MessageBox.Show("Cannot find column \"Experiments\" in the ProteinMergedResults.txt\n");
              return;
          }
          if(IfCalculateiBAQ.IsChecked==true)
          {
              if (bIfFindiBAQColumn == false)
              {
                  System.Windows.MessageBox.Show("Cannot find column \"iBAQ\" in the ProteinMergedResults.txt\n");
                  return;
              }
          }
          else
          {
              iiBAQColumnNum = -1;
          }
          if(bIfFindLFAQpepColumn == false)
          {
              System.Windows.MessageBox.Show("Cannot find column \"LFAQ\" in the ProteinMergedResults.txt\n");
              return;
          }
          if(IfCalculateTop3.IsChecked==true)
          {
              if (bIfFindTop3Column == false)
              {
                  System.Windows.MessageBox.Show("Cannot find column \"Top3\" in the ProteinMergedResults.txt\n");
                  return;
              }
          }
          else
          {
              iTop3ColumnNum = -1;
          }
          if(bIfFindPredictMolOfLFAQColumn == false)
          {
              System.Windows.MessageBox.Show("Cannot find column \"PredictedMol(LFAQ)\" in the ProteinMergedResults.txt\n");
              return;
          }
          if (IfCalculateiBAQ.IsChecked == true)
          {
              if (bIfFindPredictMolOfiBAQColumn == false)
              {
                  System.Windows.MessageBox.Show("Cannot find column \"PredictedMol(iBAQ)\" in the ProteinMergedResults.txt\n");
                  return;
              }
          }
          if(IfCalculateTop3.IsChecked==true)
          {
              if(bIfFindPredictMolOfTop3Column==false)
              {
                  System.Windows.MessageBox.Show("Cannot find column \"PredictedMol(Top3)\" in the ProteinMergedResults.txt\n");
                  return;
              }
          }
          if(bIfFindPeptidesIntensitiesColumn == false)
          {
              System.Windows.MessageBox.Show("Cannot find column \"PeptidesIntensities\" in the ProteinMergedResults.txt\n");
              return;
          }

        }
        private void GetExperimentNames(string ProteinLFAQpepPath, out Dictionary<string, int> ExperimentNamesAndIndex)
        {
            try
            {
                ExperimentNamesAndIndex = new Dictionary<string, int>();
                using (StreamReader streamReader = new StreamReader(ProteinLFAQpepPath, Encoding.Default))
                {
                    string strLine = "";
                    int iProteinIdColumnNum = 0;
                    int iMajorProteinIdColumnNum = 0;
                    int iExperimentColumnNum = 0;
                    int iiBAQColumnNum = 0;
                    int iLFAQpepColumnNum = 0;
                    int iTop3ColumnNum = 0;
                    int iPredictMolOfLFAQColumnNum = 0;
                    int iPredictMolOfiBAQColumnNum = 0;
                    int iPredictMolOfTop3ColumnNum = 0;
                    int iPeptidesIntensitiesColumnNum = 0;
                    string strExperimentNamesTemp = "";
                    string strMaxExperimentNamesTemp = "";
                    int iMaxExperimentLength = 0;
                    int iExperimentLength = 0;


                    strLine = streamReader.ReadLine();
                    string[] strValuesTemp = strLine.Split('\t');
                    int iColumnNum = strValuesTemp.Length;
                    ParseHeadLine(strValuesTemp, out  iProteinIdColumnNum, out  iMajorProteinIdColumnNum,
                                 out  iExperimentColumnNum, out  iiBAQColumnNum, out  iLFAQpepColumnNum,
                                 out  iTop3ColumnNum, out  iPredictMolOfLFAQColumnNum, out  iPredictMolOfiBAQColumnNum,
                                 out  iPredictMolOfTop3ColumnNum, out  iPeptidesIntensitiesColumnNum);

                    strLine = streamReader.ReadLine();
                    strValuesTemp = strLine.Split('\t');
                    if (strValuesTemp.Length != iColumnNum)
                    {
                        System.Windows.MessageBox.Show("The column number of line " +
                            strLine + " does not equal to the column number of the head line");
                        return;
                    }
                    while (strLine != null)
                    {
                        strValuesTemp = strLine.Split('\t');
                        if (strValuesTemp.Length != iColumnNum)
                        {
                            System.Windows.MessageBox.Show("The column number of line " +
                                strLine + " does not equal to the column number of the head line");
                            return;
                        }
                        strExperimentNamesTemp = strValuesTemp[iExperimentColumnNum];
                        iExperimentLength = strExperimentNamesTemp.Length;
                        if (iExperimentLength > iMaxExperimentLength)
                        {
                            iMaxExperimentLength = iExperimentLength;
                            strMaxExperimentNamesTemp = strExperimentNamesTemp;
                        }
                        strLine = streamReader.ReadLine();
                    }

                    string[] ExperimentNamesArrayTemp = new string[1];
                    ExperimentNamesArrayTemp = strMaxExperimentNamesTemp.Split(';');

                    int index = 0;
                    foreach (var item in ExperimentNamesArrayTemp)
                    {
                        ExperimentNamesAndIndex.Add(item, index);
                        index++;
                    }
                }

            }           
             catch(FormatException fe)
                {

                    System.Windows.MessageBox.Show("A format exception has been thrown!" + fe.Message);
                    ExperimentNamesAndIndex = new Dictionary<string, int>();
                    return;
                }

        }
        private void Global_tab_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            if (e.Source is System.Windows.Controls.TabControl)
            {
                if (VisualizationTabItem.IsSelected)
                {
                    if(IfCalculateiBAQ.IsChecked==false)
                    {
                        LFAQ_result_grid.Columns[4].Visibility = Visibility.Collapsed;
                        LFAQ_result_grid.Columns[7].Visibility = Visibility.Collapsed;
                        dataGrid_StandProteins.Columns[5].Visibility = Visibility.Collapsed;
                        dataGrid_StandProteins.Columns[8].Visibility = Visibility.Collapsed;
                    }
                    if (IfCalculateiBAQ.IsChecked == true)
                    {
                        LFAQ_result_grid.Columns[4].Visibility = Visibility.Visible;
                        LFAQ_result_grid.Columns[7].Visibility = Visibility.Visible;
                        dataGrid_StandProteins.Columns[5].Visibility = Visibility.Visible;
                        dataGrid_StandProteins.Columns[8].Visibility = Visibility.Visible;
                    }
                   if(IfCalculateTop3.IsChecked==false)
                   {
                       LFAQ_result_grid.Columns[5].Visibility = Visibility.Collapsed;
                       LFAQ_result_grid.Columns[8].Visibility = Visibility.Collapsed;
                       dataGrid_StandProteins.Columns[6].Visibility = Visibility.Collapsed;
                       dataGrid_StandProteins.Columns[9].Visibility = Visibility.Collapsed;
                   }
                   if (IfCalculateTop3.IsChecked == true)
                   {
                       LFAQ_result_grid.Columns[5].Visibility = Visibility.Visible;
                       LFAQ_result_grid.Columns[8].Visibility = Visibility.Visible;
                       dataGrid_StandProteins.Columns[6].Visibility = Visibility.Visible;
                       dataGrid_StandProteins.Columns[9].Visibility = Visibility.Visible;
                   }      
                    if (bIfQuanted)
                    {
                        Dictionary<string, double> StandProteins;
                        LoadStandProteins(out StandProteins);
                        ComboBoxItem cbiIfContainStand = (ComboBoxContainStand.SelectedItem as ComboBoxItem);
                       
                        try
                        {
                            ProteinsResult = new ObservableCollection<ProteinResult>();
                            StandProteinsResult = new ObservableCollection<StandProteinResult>();
                            SelectExperiments = new ObservableCollection<ExperimentName>();
                            string ProteinLFAQpepPath = txtResultPath.Text + "\\ProteinMergedResults.txt";
                            bool IfNewExperiment = true;
                            Dictionary<string, int> ExperimentNamesAndIndex;
                            ExperimentNamesAndIndex = new Dictionary<string, int>();
                            GetExperimentNames(ProteinLFAQpepPath,out ExperimentNamesAndIndex);

                            using (StreamReader streamReader = new StreamReader(ProteinLFAQpepPath, Encoding.Default))
                            {
                                string strLine = "";
                                string strProteinNameTemp;
                                //string strExperimentTemp;
                                double dLFAQpepTemp = 0.0;
                                double dMaxquantIBAQTemp = 0.0;
                                double dTop3Temp=0.0;
                                double dPredictedMolOfiBAQTemp = 0.0;
                                double dPredictedMolOfTop3Temp = 0.0;
                                double dPredictedMolOfLFAQpepTemp = 0.0;
                                int iExperiment = 0;
                                double dpeptideIntensity;
                                bool bIfFindMolsOfStandProtein;

                               int iProteinIdColumnNum;
                               int iMajorProteinIdColumnNum;
                               int iExperimentColumnNum;
                               int iiBAQColumnNum;
                               int iLFAQpepColumnNum; 
                               int iTop3ColumnNum;
                               int iPredictMolOfLFAQColumnNum;
                               int iPredictMolOfiBAQColumnNum;
                               int iPredictMolOfTop3ColumnNum;
                               int iPeptidesIntensitiesColumnNum;
                               //int iBegin, iFinish;
                               string[]LFAQArrayTemp=new string[1];
                               string[] iBAQArrayTemp=new string[1];
                               string[] Top3ArrayTemp=new string[1];
                               string[] PredictedMolsOfLFAQsTemp = new string[1];
                               string[] PredictedMolsOfiBAQsTemp = new string[1];
                               string[] PredictedMolsOfTop3sTemp = new string[1];
                               string[] PeptidesIntensitiesTemp = new string[1];

                               string strExperimentsTemp;
                               string[] ExperimentsTemp = new string[1];
                               ArrayList ExperimentNameIndices=new ArrayList();

                               string strTemp;

                                string strProteinIDsTemp;
                                string strPeptidesIntensities;
                                int iPeptidesBegin, iPeptidesEnd;

                                strLine = streamReader.ReadLine();
                                string []strValuesTemp = strLine.Split('\t');
                                int iColumnNum=strValuesTemp.Length;
                                ParseHeadLine(strValuesTemp, out  iProteinIdColumnNum, out  iMajorProteinIdColumnNum,
                                             out  iExperimentColumnNum, out  iiBAQColumnNum, out  iLFAQpepColumnNum, 
                                             out  iTop3ColumnNum,out  iPredictMolOfLFAQColumnNum,out  iPredictMolOfiBAQColumnNum,
                                             out  iPredictMolOfTop3ColumnNum,out  iPeptidesIntensitiesColumnNum);

                                strLine = streamReader.ReadLine();
                                strValuesTemp = strLine.Split('\t');
                                if(strValuesTemp.Length!=iColumnNum)
                                {
                                    System.Windows.MessageBox.Show("The column number of line " +
                                        strLine + " does not equal to the column number of the head line");
                                    return;
                                }
                                List<string> MissStandProteins = new List<string>();
                                MissStandProteins.Add("Warning: Cannot find these proteins in StandardProteins.txt:");
                                while (strLine != null)
                                {
                                    strValuesTemp = strLine.Split('\t');
                                    if (strValuesTemp.Length != iColumnNum)
                                    {
                                        System.Windows.MessageBox.Show("The column number of line " +
                                            strLine + " does not equal to the column number of the head line");
                                        return;
                                    }

                                    strProteinIDsTemp = strValuesTemp[iProteinIdColumnNum];
                                    strProteinNameTemp = strValuesTemp[iMajorProteinIdColumnNum];
                                    strExperimentsTemp = strValuesTemp[iExperimentColumnNum];
                                    ExperimentsTemp = strExperimentsTemp.Split(';');
                                    ExperimentNameIndices.Clear();
                                    foreach (var item in ExperimentsTemp)
                                    {
                                        ExperimentNameIndices.Add(ExperimentNamesAndIndex[item]);
                                    }
                                    iExperiment = 0;
                                    foreach (int Experimentindex in ExperimentNameIndices)
                                    {
                                        ProteinResult ProteinResultTemp = new ProteinResult();
                                        ProteinResultTemp.ListPeptidesIntensity = new List<double>();
                                        ProteinResultTemp.ListPeptidesIntensity.Clear();

                                        if (IfCalculateiBAQ.IsChecked == true)
                                        {
                                            strTemp = strValuesTemp[iiBAQColumnNum];  //iBAQ
                                            iBAQArrayTemp = strTemp.Split(';');
                                            dMaxquantIBAQTemp = Convert.ToDouble(iBAQArrayTemp[Experimentindex]);
                                            ProteinResultTemp.MaxquantIBAQ = dMaxquantIBAQTemp;

                                            strTemp = strValuesTemp[iPredictMolOfiBAQColumnNum];
                                            PredictedMolsOfiBAQsTemp = strTemp.Split(';');
                                            dPredictedMolOfiBAQTemp = Convert.ToDouble(PredictedMolsOfiBAQsTemp[Experimentindex]);
                                            ProteinResultTemp.m_PredictedMolOfiBAQ = dPredictedMolOfiBAQTemp;
                                        }

                                        strTemp = strValuesTemp[iLFAQpepColumnNum];  //LFAQpep
                                        LFAQArrayTemp = strTemp.Split(';');
                                        dLFAQpepTemp = Convert.ToDouble(LFAQArrayTemp[Experimentindex]);
                                        ProteinResultTemp.LFAQpep = dLFAQpepTemp;

                                        if (IfCalculateTop3.IsChecked == true)
                                        {
                                            strTemp = strValuesTemp[iTop3ColumnNum];
                                            Top3ArrayTemp = strTemp.Split(';');
                                            dTop3Temp = Convert.ToDouble(Top3ArrayTemp[Experimentindex]);
                                            ProteinResultTemp.Top3 = dTop3Temp;

                                            strTemp = strValuesTemp[iPredictMolOfTop3ColumnNum];
                                            PredictedMolsOfTop3sTemp = strTemp.Split(';');
                                            dPredictedMolOfTop3Temp = Convert.ToDouble(PredictedMolsOfTop3sTemp[Experimentindex]);
                                            ProteinResultTemp.m_PredictedMolOfTop3 = dPredictedMolOfTop3Temp;

                                        }

                                        strTemp = strValuesTemp[iPredictMolOfLFAQColumnNum];
                                        PredictedMolsOfLFAQsTemp = strTemp.Split(';');
                                        dPredictedMolOfLFAQpepTemp = Convert.ToDouble(PredictedMolsOfLFAQsTemp[Experimentindex]);
                                        ProteinResultTemp.PredictedMolOfLFAQpep = dPredictedMolOfLFAQpepTemp;

                                        strTemp = strValuesTemp[iPeptidesIntensitiesColumnNum];
                                        PeptidesIntensitiesTemp = strTemp.Split(';');
                                        strPeptidesIntensities = PeptidesIntensitiesTemp[Experimentindex];

                                        iPeptidesBegin = 0;
                                        iPeptidesEnd = strPeptidesIntensities.IndexOf(",", iPeptidesBegin);
                                        while (iPeptidesEnd != -1)
                                        {
                                            dpeptideIntensity = Convert.ToDouble(strPeptidesIntensities.Substring(iPeptidesBegin, iPeptidesEnd - iPeptidesBegin));
                                            ProteinResultTemp.ListPeptidesIntensity.Add(dpeptideIntensity);
                                            iPeptidesBegin = iPeptidesEnd + 1;
                                            iPeptidesEnd = strPeptidesIntensities.IndexOf(",", iPeptidesBegin);
                                        }
                                        dpeptideIntensity = Convert.ToDouble(strPeptidesIntensities.Substring(iPeptidesBegin, strPeptidesIntensities.Length - iPeptidesBegin));
                                        ProteinResultTemp.ListPeptidesIntensity.Add(dpeptideIntensity);

                                        // these code will make sure the set contains the same item
                                        ProteinResultTemp.strProteinName = strProteinNameTemp;
                                        ProteinResultTemp.strProteinIDs = strProteinIDsTemp;
                                        ProteinResultTemp.strExperiment = ExperimentsTemp[iExperiment];
                                        ProteinsResult.Add(ProteinResultTemp);
                                       
                                        if (cbiIfContainStand.Content.ToString() == "yes")
                                        {
                                            if (strProteinNameTemp.IndexOf(txtIdentifierOfStand.Text) != -1)
                                            {
                                                bIfFindMolsOfStandProtein = false;
                                                foreach (var item in StandProteins)
                                                {
                                                    if (strProteinNameTemp.IndexOf(item.Key) != -1)
                                                    {
                                                        StandProteinsResult.Add(new StandProteinResult
                                                        {
                                                            strProteinName = strProteinNameTemp,
                                                            strProteinIDs = strProteinIDsTemp,
                                                            strExperiment = ExperimentsTemp[iExperiment],
                                                            MaxquantIBAQ = dMaxquantIBAQTemp,
                                                            Top3 = dTop3Temp,
                                                            LFAQpep = dLFAQpepTemp,
                                                            PredictedMolOfLFAQpep = dPredictedMolOfLFAQpepTemp,
                                                            PredictedMolOfiBAQ = dPredictedMolOfiBAQTemp,
                                                            PredictedMolOfTop3 = dPredictedMolOfTop3Temp,
                                                            SpikedInMols = item.Value

                                                        });
                                                        bIfFindMolsOfStandProtein = true;
                                                    }
                                                    IfNewExperiment = true;
                                                    foreach (ExperimentName Exp in SelectExperiments)
                                                    {
                                                        if (Exp.strExperimentName == ExperimentsTemp[iExperiment])
                                                        {
                                                            IfNewExperiment = false;
                                                        }
                                                    }
                                                    if (IfNewExperiment)
                                                    {

                                                        SelectExperiments.Add(new ExperimentName
                                                        {
                                                            strExperimentName = ExperimentsTemp[iExperiment]
                                                        });
                                                    }

                                                }
                                            }
                                        }
                                        iExperiment++;
                                    }// end foreach ExperimentNameIndices)
                                    strLine = streamReader.ReadLine();
                                }
  
                                LFAQ_result_grid.DataContext = ProteinsResult;
                                if (cbiIfContainStand.Content.ToString() == "yes")
                                {
                                    dataGrid_StandProteins.DataContext = StandProteinsResult;
                                }
                                else
                                {
                                    dataGrid_StandProteins.DataContext = "";
                                }
                                DoEvents();
                                if(MissStandProteins.Count>1)
                                {
                                    var message = String.Join(Environment.NewLine, MissStandProteins);
                                    System.Windows.MessageBox.Show(message);
                                }
                                streamReader.Close();
                            }  // end using stream

                        }// end try
                        catch (IOException ex)
                        {
                            LFAQ_result_grid.DataContext = ProteinsResult;
                            if (cbiIfContainStand.Content.ToString() == "yes")
                            {
                                dataGrid_StandProteins.DataContext = StandProteinsResult;
                                dataGrid_SelectExperiment.DataContext = SelectExperiments;
                            }
                            LoadEmptyLinearFittingPic();
                            LoadEmptyBarPicOfOneProtein();
                            System.Windows.MessageBox.Show("An IOException has been thrown!" + ex.Message);
                            return;
                        }
                        catch(FormatException fe)
                        {

                            System.Windows.MessageBox.Show("A format exception has been thrown!" + fe.Message);
                            return;
                        }
                        if (cbiIfContainStand.Content.ToString() == "yes")
                        {
                            
                            string StandProteinsResultsPath ="StandardProteinResults.txt";
                            SaveStandProteinsResults(StandProteinsResult, StandProteinsResultsPath);
                            bIfSaveStandProteins = true;
                            string CurrentDirectory = System.Environment.CurrentDirectory;
                            CreatLinearFittingScript(CurrentDirectory, StandProteinsResultsPath);
                            dataGrid_SelectExperiment.DataContext = SelectExperiments;
                        }
                        else
                        {
                            LoadEmptyLinearFittingPic();
                        }
                    }// end if (bIfQuanted)
                    else
                    {
                        LFAQ_result_grid.DataContext = null;
                        LoadEmptyBarPicOfOneProtein();
                        dataGrid_StandProteins.DataContext = null;
                        LoadEmptyLinearFittingPic();
                    }
                }// end if(VisualizationTabItem.IsSelected)
            }

        }


        private void LoadStandProteins( out Dictionary<string, double> StandProteins)
        {
            StandProteins = new Dictionary<string, double>();
            StandProteins.Clear();
            string strLine;
            int iStart, iEnd;
            string CurrentDirectory = System.Environment.CurrentDirectory;
            string StandProteinsPath = txtResultPath.Text+"\\StandardProteins.txt";
            if(!File.Exists(StandProteinsPath))
            {
                StandProteinsPath = "StandardProteins.txt";
            }
            string strProteinName;
            double dSpikedInOfProtein;
            //Dictionary<string, double> StandProteins = new Dictionary<string, double>();
            using (StreamReader streamReader = new StreamReader(StandProteinsPath, Encoding.Default))
            {
                strLine = streamReader.ReadLine();
                strLine = streamReader.ReadLine();
                while (strLine != null)
                {
                    iStart = 0;
                    iEnd = strLine.IndexOf("\t");
                    strProteinName = strLine.Substring(iStart, iEnd - iStart);  //uniProt Accession Number

                    //iStart = iEnd + 1;
                    //iEnd = strLine.IndexOf("\t", iStart);  //UPS1 Amount (fmol)
                    //if (iEnd == -1)
                    //{
                    //    iEnd = strLine.IndexOf("\n", iStart);    
                    //}
                    iStart = iEnd + 1;
                    iEnd = strLine.IndexOf("\t", iStart);   //UPS2 Amount (fmol)
                    if (iEnd == -1)
                    {
                        iEnd = strLine.IndexOf("\n", iStart);
                    }
                    //////###
                    if(iEnd != -1)
                    {
                        dSpikedInOfProtein = Convert.ToDouble(strLine.Substring(iStart, iEnd - iStart));  
                    }
                    else
                    {
                        dSpikedInOfProtein = Convert.ToDouble(strLine.Substring(iStart, strLine.Length - iStart));
                    }

                    

                    try
                    {
                        if(StandProteins.ContainsKey(strProteinName)==false)
                        {
                            StandProteins.Add(strProteinName, dSpikedInOfProtein);
                        }

                    }
                    catch (System.Exception excep)
                    {
                        System.Windows.MessageBox.Show(excep.ToString());
                    }

                    strLine = streamReader.ReadLine();
                }
            }
         }
        public void SaveStandProteinsResults(ObservableCollection<StandProteinResult> standProteinsResults, string path)
        {
            try
            {          
                StreamWriter writer = File.CreateText(path);
                writer.WriteLine("ProteinName\tExperiemntName\tSpikedInMols\tLFAQ\t");
                foreach(var item in standProteinsResults)
                {
                    if (item.LFAQpep == 0.0)
                    {
                        writer.WriteLine(item.strProteinName + "\t"+item.strExperiment+"\t" + Math.Log10(item.SpikedInMols) + "\t0.0\t");
                    }
                    else
                    {
                        writer.WriteLine(item.strProteinName + "\t" + item.strExperiment + "\t" + Math.Log10(item.SpikedInMols) + "\t" + Math.Log10(item.LFAQpep) + "\t");
                    }
                }
                writer.Close();
            }
            catch (IOException ex)
            {
                System.Windows.MessageBox.Show("An IOException has been thrown!"+ex.Message);
                return;
            }
        }
        public void CreatLinearFittingScript(string workPath=@"D:",string dataPath="" )
        {
            string scriptPath = System.IO.Path.GetFullPath(workPath) + @"\LinearFitOfLFAQ.r";
            workPath = workPath.Replace("\\", "\\\\");

            if (dataGrid_SelectExperiment.SelectedIndex == -1)
                dataGrid_SelectExperiment.SelectedIndex = 0;
            if(dataGrid_SelectExperiment.SelectedIndex>=0&&dataGrid_SelectExperiment.SelectedIndex<SelectExperiments.Count)
            {
                ExperimentName experimentName = SelectExperiments[dataGrid_SelectExperiment.SelectedIndex];
                dataPath = dataPath.Replace("\\", "\\\\");
                string str = "";
                using (StreamWriter strW = new StreamWriter(scriptPath, false, Encoding.GetEncoding("GB2312")))
                {
                    str = "setwd(\"" + workPath + "\")\n";
                    str += "library(ggplot2)\n";
                    str += "d<-read.table(\"" + dataPath + "\",sep=\"\\t\",header = TRUE)\n";
                    str += "x<-d[d$ExperiemntName==\"" + experimentName.strExperimentName + "\",\"SpikedInMols\"];\n";
                    str += "LFAQ<-d[d$ExperiemntName==\"" + experimentName.strExperimentName + "\",\"LFAQ\"];\n";
                   // str += "LFAQpro<-d[d$ExperiemntName==\"" + experimentName.strExperimentName + "\",\"LFAQpro\"];\n";
                    str += "LLFAQ<-rep(\"LFAQ\",length(LFAQ))\n";
                    //str += "LLFAQpro<-rep(\"LFAQpro\",length(LFAQpro))\n";
                    //str += "plotFrame<-data.frame(SpikedIn=c(x,x),Label=c(LLFAQpep,LLFAQpro),LFAQ=c(LFAQpep,LFAQpro))\n";
                    str += "plotFrame<-data.frame(SpikedIn=x,Label=LLFAQ,LFAQ=LFAQ)\n";
                    str += "tiff(file=\"" + workPath + "/LFAQOfStandPro.tiff\",res=500,width=2000,height=1500,compression = \"lzw\")\n";
                    str += "ggplot(plotFrame,aes(SpikedIn,LFAQ,colour=Label))+geom_point()+geom_smooth(method = \"lm\")+ylab(\"LFAQ protein intensity(log10)\")+xlab(\"Spiked-in standards amount(log10)\")\n";
                    str += "dev.off()\n";
                    strW.WriteLine(str);
                    strW.Close();
                }
            }
            

            
            
        }
        private void LinearFittingRun(string CurrentDirectory)
        {
            //Process LoadProcess = new Process();
            string ExeName = "Rscript.exe";
            string args = "\""+CurrentDirectory + "\\LinearFitOfLFAQ.r\"";
            ProcessStartInfo StartInfo = new ProcessStartInfo(ExeName, args);
            StartInfo.UseShellExecute = false;
            StartInfo.CreateNoWindow = true;
            //hide  the window of exe
            //StartInfo.WindowStyle = ProcessWindowStyle.Hidden;
            try
            {
                // run exe
                Process proBach = Process.Start(StartInfo);
                while (!proBach.HasExited)
                {
                    proBach.WaitForExit();
                }
            }
            catch(Win32Exception e)
            {
                System.Windows.MessageBox.Show("Please check if R has been installed correctly, see user guide of LFAQ for details.");
                bIfInstallR = false;
            }
            catch (InvalidOperationException ive)
            {
                System.Windows.MessageBox.Show("There are something wrong with the R program \"LinearFitOfLFAQ.r\".");
                bIfInstallR = false;
                return;
            }

        }
        private void LoadLinearFittingPic()
        {
          try
          {
              BitmapImage bi = new BitmapImage();
              // BitmapImage.UriSource must be in a BeginInit/EndInit block.
              bi.BeginInit();

              FileStream fs = new FileStream("LFAQOfStandPro.tiff", FileMode.Open);
              byte[] byData = new byte[fs.Length];
              fs.Read(byData, 0, byData.Length);
              fs.Close();
              bi.StreamSource = new MemoryStream(byData);

              bi.EndInit();
              // Set the image source.
              LFAQ_image.Source = bi;
          }
              catch (ArgumentException)
          {
              System.Windows.MessageBox.Show("Please make sure the R has been installed successfully.");
          }
              catch(EndOfStreamException)
          {
              System.Windows.MessageBox.Show("The LFAQOfStandPro.tiff is empty!");

          }
            catch (IOException ex)
            {
                //System.Windows.MessageBox.Show(ex.Message);
            }

        }
        private void LoadEmptyLinearFittingPic()
        {
            try
            {
                LFAQ_image.Source = null;
            }
            catch (IOException ex)
            {
                System.Windows.MessageBox.Show(ex.Message);
            }
        }
        private void LFAQpep_gridSelectedChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            // for OneProtein_image
            string PeptidesOfOneProteinPath = "PeptidesOfOneProtein.txt";
            if (LFAQ_result_grid.SelectedIndex == -1)
                LFAQ_result_grid.SelectedIndex = 0;
            if(LFAQ_result_grid.SelectedIndex>=0&&LFAQ_result_grid.SelectedIndex<ProteinsResult.Count)
            {
                ProteinResult OneProtein = ProteinsResult[LFAQ_result_grid.SelectedIndex];
                SavepeptidesOfOneProtein(OneProtein, PeptidesOfOneProteinPath);
                string CurrentDirectory = System.Environment.CurrentDirectory;
                CreatPicOfOneProteinScript(CurrentDirectory, PeptidesOfOneProteinPath);
                if(bIfInstallR)
                {
                    PicOfOneProteinRun(CurrentDirectory);
                    LoadBarPicOfOneProtein();
                }

            }
           

        }
      private void  SavepeptidesOfOneProtein(ProteinResult OneProtein,string path)
        {
            try
            {
                StreamWriter writer = File.CreateText(path);
                writer.WriteLine("PeptideIndex\tPeptideIntensity\t");
                int i=1;
                foreach (var item in OneProtein.ListPeptidesIntensity)
                {
                    writer.WriteLine(i + "\t" + item +"\t");
                    i++;
                }
                writer.Close();
            }
            catch (IOException ex)
            {
                System.Windows.MessageBox.Show("An IOException has been thrown!"+ex.Message);
                return;
            }

        }
      private void CreatPicOfOneProteinScript(string workPath=@"D:",string dataPath="")
      {
          string scriptPath = System.IO.Path.GetFullPath(workPath) + @"\BarOfOneProtein.r";
          workPath = workPath.Replace("\\", "\\\\");

          dataPath = dataPath.Replace("\\", "\\\\");
          string str = "";
          using (StreamWriter strW = new StreamWriter(scriptPath, false, Encoding.GetEncoding("GB2312")))
          {
              str+= "a<-c(1,2)\n";
              str+= "write.csv(a,\'Intensities0.csv\')\n";
              str = "setwd(\"" + workPath + "\")\n";
              str += "library(ggplot2)\n";
              str += "d<-read.table(\"" + dataPath + "\",sep=\"\\t\",header = TRUE)\n";
              str += "Peptidesindex<-d[,1];\n";
              str += "Intensities<-d[,2];\n";
              //str += "write.csv(Intensities,\'Intensities.csv\');\n";
              str += "plotFrame<-data.frame(Peptidesindex,Intensities);\n";
              str += "tiff(file=\"" + workPath + "/BarOfOneProtein.tiff\",res=500,width=2000,height=1500,compression = \"lzw\")\n";
              str += "if(length(Peptidesindex)<10){\n";
              str += "ggplot()+geom_bar(aes(Peptidesindex,Intensities),stat = \"identity\",fill=Peptidesindex)+ ylab(\"LFAQ peptide intensity\")+ xlab(\"Peptide index\")+scale_x_continuous(breaks=Peptidesindex)\n";
              str += "}else{\n";
              str += "ggplot()+geom_bar(aes(Peptidesindex,Intensities),stat = \"identity\",fill=Peptidesindex)+ ylab(\"LFAQ peptide intensity\")+ xlab(\"Peptide index\")\n";
              str += "}\n";
              str += "dev.off()\n";

              strW.WriteLine(str);
              strW.Close();
          }
      }

      private void PicOfOneProteinRun(string CurrentDirectory)
      {
          //Process LoadProcess = new Process();
          string ExeName = "Rscript.exe";

          string args = "\""+CurrentDirectory + "\\BarOfOneProtein.r\"";
          ProcessStartInfo StartInfo = new ProcessStartInfo(ExeName, args);
          StartInfo.UseShellExecute = false;
          StartInfo.CreateNoWindow = true;
          //hide  the window of exe
         // StartInfo.WindowStyle = ProcessWindowStyle.Hidden;
          // run exe
          try 
          {
              Process proBach = Process.Start(StartInfo);
              while (!proBach.HasExited)
              {
                  proBach.WaitForExit();
              }
              if(proBach.ExitCode!=0)
              {
                  System.Windows.MessageBox.Show("There are something wrong with the R script \"BarOfOneProtein.r\".");
                  return;
              }
          }
          catch (InvalidOperationException ive)
          {
              System.Windows.MessageBox.Show("There are something wrong with the R program \"BarOfOneProtein.r\".");
              bIfInstallR = false;
              return;
          }
          catch (Win32Exception ex)
          {
              System.Windows.MessageBox.Show("Please check if R has been installed correctly, see user guide of LFAQ for details.");
              bIfInstallR = false;
              return;
          }

      }
      private void LoadBarPicOfOneProtein()
      {
          try
          {
              BitmapImage bi = new BitmapImage();
              // BitmapImage.UriSource must be in a BeginInit/EndInit block.
              bi.BeginInit();

              FileStream fs = new FileStream("BarOfOneProtein.tiff", FileMode.Open);
              byte[] byData = new byte[fs.Length];
              fs.Read(byData, 0, byData.Length);
              fs.Close();
              bi.StreamSource = new MemoryStream(byData);

              bi.EndInit();
              // Set the image source.
              OneProtein_image.Source = bi;

          }
          catch (ArgumentException)
          {
              System.Windows.MessageBox.Show("BarOfOneProtein.tiff is empty!");
          }
              catch(EndOfStreamException)
          {
              System.Windows.MessageBox.Show("Please make sure the R has been installed successfully.");
          }
          catch (IOException ex)
          {
              //System.Windows.MessageBox.Show(ex.ToString());
              //System.Windows.MessageBox.Show(ex.Message);
              System.Windows.MessageBox.Show("Cannot open BarOfOneProtein.tiff");
              return;
          }

      }
      private void LoadEmptyBarPicOfOneProtein()
      {
          try
          {
            // Set the image source.
              OneProtein_image.Source = null;
          }
          catch (IOException ex)
          {
              System.Windows.MessageBox.Show(ex.Message);
              return;
          }

      }
        private BitmapImage BitmapToBitmapImage(System.Drawing.Bitmap bitmap)
        {
            BitmapImage bitmapImage = new BitmapImage();

            using (System.IO.MemoryStream ms = new System.IO.MemoryStream())
            {
                ms.Position = 0;
                //bitmap.Save(ms, bitmap.RawFormat);
                bitmap.Save(ms, System.Drawing.Imaging.ImageFormat.Png);
                bitmapImage.BeginInit();
                bitmapImage.StreamSource = ms;
                bitmapImage.CacheOption = BitmapCacheOption.OnLoad;
                bitmapImage.EndInit();
                bitmapImage.Freeze();
            }

            return bitmapImage;
        }

        private void ShowError()
        {
           string logPath = txtResultPath.Text + "\\log.txt";
            if(!File.Exists(logPath))
            {
                System.Windows.MessageBox.Show("Cannot find file: " + logPath);
                return;
            }
            try
            {
                using (StreamReader sr = new StreamReader(logPath, Encoding.Default))
                {
                    string Line;
                    string strTemp;
                    int iEnd;
                    while ((Line = sr.ReadLine()) != null)
                    {
                        iEnd = Line.IndexOf("\t");
                        if (iEnd != -1)
                        {
                            strTemp = Line.Substring(0, iEnd);
                            if (strTemp == "Error:")
                            {
                                System.Windows.MessageBox.Show(Line);
                                break;
                            }
                        }
                    }
                }

            }
            catch (IOException ex)
            {
                System.Windows.MessageBox.Show("An IOException has been thrown!"+ex.Message);
                return;
            }
        }

        public class ProteinResult
        {
            public event PropertyChangedEventHandler PropertyChanged;
            private string m_strProteinName;
            private string m_strProteinIDs;
            private string m_experimentName;
            private double m_MaxquantIBAQ;
            private double m_dTop3;
            private double m_LFAQpep;
            public double m_PredictedMolOfLFAQpep;
            public double m_PredictedMolOfiBAQ;
            public double m_PredictedMolOfTop3;
            public string strExperiment
            {
                get { return m_experimentName; }
                set
                {
                    m_experimentName = value;
                    if (PropertyChanged != null)
                    {
                        PropertyChanged(this, new PropertyChangedEventArgs("strExperiment"));
                    }
                }
            }

            public string strProteinName
            {
                get { return m_strProteinName; }
                set
                {
                    m_strProteinName = value;
                    if (PropertyChanged != null)
                    {
                        PropertyChanged(this, new PropertyChangedEventArgs("strProteinName"));
                    }
                }
            }

            public string strProteinIDs
            {
                get { return m_strProteinIDs; }
                set
                {
                    m_strProteinIDs = value;
                    if (PropertyChanged != null)
                    {
                        PropertyChanged(this, new PropertyChangedEventArgs("strProteinIDs"));
                    }
                }
            }

            public double MaxquantIBAQ
            {
                get { return m_MaxquantIBAQ; }
                set
                {
                    m_MaxquantIBAQ = value;
                    if (PropertyChanged != null)
                    {
                        PropertyChanged(this, new PropertyChangedEventArgs("MaxQuantIBAQ"));
                    }
                }
            }
            public double Top3
            {
                get { return m_dTop3; }
                set
                {
                    m_dTop3 = value;
                    if (PropertyChanged != null)
                    {
                        PropertyChanged(this, new PropertyChangedEventArgs("Top3"));
                    }
                }
            }
            public double LFAQpep
            {
                get { return m_LFAQpep; }
                set
                {
                    m_LFAQpep = value;
                    if (PropertyChanged != null)
                    {
                        PropertyChanged(this, new PropertyChangedEventArgs("LFAQpep"));
                    }
                }
            }

           public double PredictedMolOfLFAQpep
            {
                get { return m_PredictedMolOfLFAQpep; }
                set
                {
                    m_PredictedMolOfLFAQpep = value;
                    if (PropertyChanged != null)
                    {
                        PropertyChanged(this, new PropertyChangedEventArgs("PredictedMolOfLFAQpep"));
                    }
                }
            }
            public double PredictedMolOfiBAQ
            {
                get { return m_PredictedMolOfiBAQ; }
                set
                {
                    m_PredictedMolOfiBAQ = value;
                    if (PropertyChanged != null)
                    {
                        PropertyChanged(this, new PropertyChangedEventArgs("PredictedMolOfiBAQ"));
                    }
                }
            }
            public double PredictedMolOfTop3
            {
                get { return m_PredictedMolOfTop3; }
                set
                {
                    m_PredictedMolOfTop3 = value;
                    if (PropertyChanged != null)
                    {
                        PropertyChanged(this, new PropertyChangedEventArgs("PredictedMolOfTop3"));
                    }
                }
            }
           
            public List<double> ListPeptidesIntensity;
        }


        public class StandProteinResult
        {
            public event PropertyChangedEventHandler PropertyChanged;
            private string m_strProteinName;
            private string m_strProteinIDs;
            private string m_experimentName;
            private double m_MaxquantIBAQ;
            private double m_dTop3;
            private double m_LFAQpep;
            //private double m_LFAQpro;
            private double m_SpikedInMols;
            private double m_PredictedMolOfLFAQpep;
           // private double m_PredictedMolOfLFAQpro;
            private double m_PredictedMolOfiBAQ;
            private double m_PredictedMolOfTop3;
            public string strExperiment
            {
                get { return m_experimentName; }
                set
                {
                    m_experimentName = value;
                    if (PropertyChanged != null)
                    {
                        PropertyChanged(this, new PropertyChangedEventArgs("strExperiment"));
                    }
                }
            }

            public string strProteinName
            {
                get { return m_strProteinName; }
                set
                {
                    m_strProteinName = value;
                    if (PropertyChanged != null)
                    {
                        PropertyChanged(this, new PropertyChangedEventArgs("strProteinName"));
                    }
                }
            }
            public string strProteinIDs
            {
                get { return m_strProteinIDs; }
                set
                {
                    m_strProteinIDs = value;
                    if (PropertyChanged != null)
                    {
                        PropertyChanged(this, new PropertyChangedEventArgs("strProteinIDs"));
                    }
                }
            }

            public double MaxquantIBAQ
            {
                get { return m_MaxquantIBAQ; }
                set
                {
                    m_MaxquantIBAQ = value;
                    if (PropertyChanged != null)
                    {
                        PropertyChanged(this, new PropertyChangedEventArgs("MaxQuantIBAQ"));
                    }
                }
            }
            public double Top3
            {
                get { return m_dTop3; }
                set
                {
                    m_dTop3 = value;
                    if (PropertyChanged != null)
                    {
                        PropertyChanged(this, new PropertyChangedEventArgs("Top3"));
                    }
                }
            }
            public double LFAQpep
            {
                get { return m_LFAQpep; }
                set
                {
                    m_LFAQpep = value;
                    if (PropertyChanged != null)
                    {
                        PropertyChanged(this, new PropertyChangedEventArgs("LFAQpep"));
                    }
                }
            }

            //public double LFAQpro
            //{
            //    get { return m_LFAQpro; }
            //    set
            //    {
            //        m_LFAQpro = value;
            //        if (PropertyChanged != null)
            //        {
            //            PropertyChanged(this, new PropertyChangedEventArgs("LFAQpro"));
            //        }
            //    }
            //}

            public double PredictedMolOfLFAQpep
            {
                get { return m_PredictedMolOfLFAQpep; }
                set
                {
                    m_PredictedMolOfLFAQpep = value;
                    if (PropertyChanged != null)
                    {
                        PropertyChanged(this, new PropertyChangedEventArgs("PredictedMolOfLFAQpep"));
                    }
                }
            }
            //public double PredictedMolOfLFAQpro
            //{
            //    get { return m_PredictedMolOfLFAQpro; }
            //    set
            //    {
            //        m_PredictedMolOfLFAQpro = value;
            //        if (PropertyChanged != null)
            //        {
            //            PropertyChanged(this, new PropertyChangedEventArgs("PredictedMolOfLFAQpro"));
            //        }
            //    }
            //}
            public double PredictedMolOfiBAQ
            {
                get { return m_PredictedMolOfiBAQ; }
                set
                {
                    m_PredictedMolOfiBAQ = value;
                    if (PropertyChanged != null)
                    {
                        PropertyChanged(this, new PropertyChangedEventArgs("PredictedMolOfiBAQ"));
                    }
                }
            }
            public double PredictedMolOfTop3
            {
                get { return m_PredictedMolOfTop3; }
                set
                {
                    m_PredictedMolOfTop3 = value;
                    if (PropertyChanged != null)
                    {
                        PropertyChanged(this, new PropertyChangedEventArgs("PredictedMolOfTop3"));
                    }
                }
            }
            public double SpikedInMols
            {
                get { return m_SpikedInMols; }
                set
                {
                    m_SpikedInMols = value;
                    if (PropertyChanged != null)
                    {
                        PropertyChanged(this, new PropertyChangedEventArgs("SpikedInMols"));
                    }
                }
            }
            

        }

        public class ExperimentName
        {
            public event PropertyChangedEventHandler PropertyChanged;
            private string m_strExperimentName;

            public string strExperimentName
            {
                get { return m_strExperimentName; }
                set
                {
                    m_strExperimentName = value;
                    if (PropertyChanged != null)
                    {
                        PropertyChanged(this, new PropertyChangedEventArgs("strExperimentName"));
                    }
                }
            }
        }
        private void mzQuantMLChecked(object sender, RoutedEventArgs e)
        {
            try
            {

                MaxQuantCheckBox.IsChecked = false;
                SWATHCheckBox.IsChecked = false;
            }
            catch (IOException ex)
            {
                System.Windows.MessageBox.Show("An IOException has been thrown!" + ex.ToString());
                return;
            }
            if (this.IsLoaded)
            {
                InputDirectoryStack.Visibility = Visibility.Hidden;
                InputFileStack.Visibility = Visibility.Visible;
            }
        }

        private void SWATHChecked(object sender, RoutedEventArgs e)
        {
            try
            {

                MaxQuantCheckBox.IsChecked = false;
                mzQuantMLCheckBox.IsChecked = false;
            }
            catch (IOException ex)
            {
                System.Windows.MessageBox.Show("An IOException has been thrown!" + ex.ToString());
                return;
            }
            if (this.IsLoaded)
            {
                InputDirectoryStack.Visibility = Visibility.Hidden;
                InputFileStack.Visibility = Visibility.Visible;
            }
        }
        private void mzQuantMLUnchecked(object sender, RoutedEventArgs e)
        {
            try
            {
                //MaxQuantCheckBox.IsChecked = true;
                //SWATHCheckBox.IsChecked = true;
            }
            catch (IOException ex)
            {
                System.Windows.MessageBox.Show("An IOException has been thrown!" + ex.ToString());
                return;
            }
        }
        private void SWATHUnchecked(object sender, RoutedEventArgs e)
        {
            try
            {
                //MaxQuantCheckBox.IsChecked = true;
                //SWATHCheckBox.IsChecked = true;
            }
            catch (IOException ex)
            {
                System.Windows.MessageBox.Show("An IOException has been thrown!" + ex.ToString());
                return;
            }
        }
        
        private void MaxQuantChecked(object sender, RoutedEventArgs e)
        {
            try
            {
                if (this.IsLoaded)
                {
                    mzQuantMLCheckBox.IsChecked = false;
                    SWATHCheckBox.IsChecked = false;
                }
            }
            catch (IOException ex)
            {
                System.Windows.MessageBox.Show("An IOException has been thrown!" + ex.ToString());
                return;
            }
            catch (NullReferenceException ex)
            {
                System.Windows.MessageBox.Show("An IOException has been thrown!" + ex.ToString());
                return;
            }
            if (this.IsLoaded)
            {
                InputDirectoryStack.Visibility = Visibility.Visible;
                InputFileStack.Visibility = Visibility.Hidden;
            }

        }

        private void MaxQuantUnchecked(object sender, RoutedEventArgs e)
        {
            try
            {
                //mzQuantMLCheckBox.IsChecked = true;
            }
            catch (IOException ex)
            {
                System.Windows.MessageBox.Show("An IOException has been thrown!" + ex.ToString());
                return;
            }
        }

        private void ContainStand_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            ComboBoxItem cbi = ((sender as System.Windows.Controls.ComboBox).SelectedItem as ComboBoxItem);

            if (cbi.Content.ToString() == "yes")
            {
                if (StandProteinGrid != null)
                {
                    StandProteinGrid.Visibility = Visibility.Visible;
                }

            }
            else if (cbi.Content.ToString() == "no")
            {
                if (StandProteinGrid != null)
                {
                    StandProteinGrid.Visibility = Visibility.Hidden;
                }

            }
        }


        private void dataGrid_SelectExperiment_SelectedChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            // for LFAQ_image
            if(bIfSaveStandProteins)
            {
                string CurrentDirectory = System.Environment.CurrentDirectory;
                string StandProteinsResultsPath = "StandardProteinResults.txt";
                CreatLinearFittingScript(CurrentDirectory, StandProteinsResultsPath);
                if(bIfInstallR)
                {
                    LinearFittingRun(CurrentDirectory);
                    LoadLinearFittingPic();
                }
            }          
        }
        public string GettxtInputPath()
        {
            if (MaxQuantCheckBox.IsChecked == true)
                return txtInputDirectoryPath.Text;
            else if (SWATHCheckBox.IsChecked == true)
            {
                return txtInputFilePath.Text;
            }
            else if (mzQuantMLCheckBox.IsChecked == true)
                return txtInputFilePath.Text.Substring(0, txtInputFilePath.Text.LastIndexOf("\\"));
            else return "null";
        }
        public string GetResultPath()
        {
            return txtResultPath.Text;
        }
        public string GetInputType()
        {
            if(MaxQuantCheckBox.IsChecked==true)
            {
                return "maxquant";
            }
            else if(mzQuantMLCheckBox.IsChecked==true)
            {
                return "mzQuantML";
            }
            else if (SWATHCheckBox.IsChecked == true)
            {
                return "SWATH";
            }
            else
            {
                return "null";
            }
        }
        public string GetIdentifierOfStand()
        {
            return txtIdentifierOfStand.Text;
        }
        private void GetTxtResultPath(out string resultPath)
        {
            resultPath=txtResultPath.Text;
        }
        private void InputFile_Browse_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openDlg = new Microsoft.Win32.OpenFileDialog();
            if(mzQuantMLCheckBox.IsChecked==true)
            {
                openDlg.Filter = "mzQuantML file| *.mzq";
            }
            else if(SWATHCheckBox.IsChecked==true)
            {
                openDlg.Filter = "SWATH file| *.csv";
            }
            if (true == openDlg.ShowDialog())
            {
                string DatafromFile = openDlg.FileName;
                txtInputFilePath.Text = DatafromFile;
            }
        }
        public void DoEvents()
        {
            DispatcherFrame frame = new DispatcherFrame();
            Dispatcher.CurrentDispatcher.BeginInvoke(DispatcherPriority.Background,
                new DispatcherOperationCallback(delegate(object f)
                {
                    ((DispatcherFrame)f).Continue = false;

                    return null;
                }
                    ), frame);
            Dispatcher.PushFrame(frame);
        }
        private void DisableParameters()
        {
            MenuFile.IsEnabled = false;
            stackPanel0.IsEnabled = false;
            stackPanel1.IsEnabled = false;
            stackPanel2.IsEnabled = false;
            GridPanel3.IsEnabled = false;
            InputDirectoryStack.IsEnabled = false;
            InputFileStack.IsEnabled = false;

            Run_buttonInData.IsEnabled = false;
            Run_ButtonInQuantification.IsEnabled = false;

            RegressionParamGroupBox.IsEnabled = false;
            QuantParamGroupBox.IsEnabled = false;
        }
        private void EnableParameters()
        {
            MenuFile.IsEnabled = true;
            stackPanel0.IsEnabled = true;
            stackPanel1.IsEnabled = true;
            stackPanel2.IsEnabled = true;
            GridPanel3.IsEnabled = true;
            InputDirectoryStack.IsEnabled = true;
            InputFileStack.IsEnabled = true;

            Run_buttonInData.IsEnabled = true;
            Run_ButtonInQuantification.IsEnabled = true;

            RegressionParamGroupBox.IsEnabled = true;
            QuantParamGroupBox.IsEnabled = true;
        }

        private void About_Click(object sender, RoutedEventArgs e)
        {
            string strAboutLFAQ = "LFAQ Version: v1.0.0\nLFAQ is an efficient tool for label-free absolute protein quantification.";
            System.Windows.MessageBox.Show(strAboutLFAQ);
        }

        private void UserGuide_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                System.Diagnostics.Process.Start("http://LFAQ.github.io/LFAQ/");
            }
            catch 
            {
            }
        }



    }

}
