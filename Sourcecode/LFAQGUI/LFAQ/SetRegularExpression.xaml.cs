using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.IO;
using System.Text.RegularExpressions;
using LFAQ;

namespace LFAQ
{
    /// <summary>
    ///Interaction logic for SetRegularExpression.xaml 
    /// </summary>

    public delegate void ChangeTextHandler(string text);
    public delegate string GetTextBoxContentHandler();
    public delegate string GetInputTypeHandle();
    public partial class SetRegularExpression : Window
    {        
        public GetTextBoxContentHandler getInputPathText;
        public GetInputTypeHandle getInputType;
        public ObservableCollection<RegularExpressionExample> exams { get; set; }
        public ObservableCollection<ProteinName> ProteinNames { get; set; }
        
        public event ChangeTextHandler ChangeTextEvent;
        public SetRegularExpression()
        {
            InitializeComponent();
            try
            {
                exams = new ObservableCollection<RegularExpressionExample>();
                exams.Add(new RegularExpressionExample { header = ">IPI:IPI00107908.1|TREMBL:Q9D5E3|ENSEMBL:ENSMUSP00000106228|REFSEQ:NP_080412", re = ">(.*?)\\|", Identifier = "IPI:IPI00107908.1" });
                exams.Add(new RegularExpressionExample { header = ">sp|P02768ups|ALBU_HUMAN_UPS Serum albumin (Chain 26-609) - Homo sapiens (Human)", re = ">sp\\|(.*?)\\|", Identifier = "P02768ups" });
                exams.Add(new RegularExpressionExample { header = ">ENSP00000381386 pep:known chromosome:GRCh37:22:24313554:24316773:-1", re = ">(.*?)\\s", Identifier = "ENSP00000381386" });
                examples.DataContext = exams;
            }
            catch (IOException ex)
            { 
                System.Windows.MessageBox.Show("An IOException has been thrown!"+ex.Message);
                return;
            }
        }

        private void fasta_browse_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openDlg = new Microsoft.Win32.OpenFileDialog();
            openDlg.Filter = "fasta Files| *.fasta";
            string strLine;
            string strProteinIdentifier;
            List<string> ProteinsInFasta = new List<string>();

            if (true == openDlg.ShowDialog()) //get protein names from fasta file
            {
                string FastaFilePath = openDlg.FileName;
                txtfasta.Text = FastaFilePath;
                int iProteins = 0;
                int iFastaProteinsNum = 999;
                try
                {
                    FileStream aFile = new FileStream(txtfasta.Text, FileMode.Open);
                    StreamReader streamReader = new StreamReader(aFile);
                    strLine = streamReader.ReadLine();
                    while (strLine != null && iProteins <= iFastaProteinsNum)
                    {
                        while (strLine == "")
                        {
                            strLine = streamReader.ReadLine();
                        }
                        if (strLine[0] == '>')
                        {
                            ProteinsInFasta.Add(strLine);
                            iProteins++;
                        }
                        strLine = streamReader.ReadLine();
                    }
                    fasta_status.Content ="Load protein headers from the fasta file.";
                    streamReader.Close();
                }
                catch (ArgumentException)
                {
                    System.Windows.MessageBox.Show("Cannot open "+txtfasta.Text);
                }
                catch (IOException ex)
                {
                    System.Windows.MessageBox.Show("An IOException has been thrown!" + ex.Message);
                    return;
                }
                int iShowLines = 100;
                string InputType = getInputType.Invoke();
                if (InputType == "maxquant")
                {
                    try
                    {
                        string strResultPath = getInputPathText.Invoke();                       
                        strResultPath = strResultPath + "\\proteinGroups.txt";
                        FileStream aFile = new FileStream(strResultPath, FileMode.Open);
                        StreamReader streamReader = new StreamReader(aFile);
                        ProteinNames = new ObservableCollection<ProteinName>();

                        iProteins = 0;
                        strLine = streamReader.ReadLine();
                        strLine = streamReader.ReadLine();
                        int iBegain, iEnd;
                        string strProteinIdTemp;
                        string strProteinLastIdTemp = "";
                        while (strLine != null&&iProteins< iShowLines)
                        {
                            while (strLine == "")
                            {
                                strLine = streamReader.ReadLine();
                            }

                            iBegain = strLine.IndexOf("\t");
                            if (iBegain == -1)
                            {
                                System.Windows.MessageBox.Show("Cannot read protien ids from proteinGroups.txt");
                            }
                            iEnd = strLine.IndexOf("\t", iBegain + 1);
                            if (iEnd == -1)
                            {
                                System.Windows.MessageBox.Show("Cannot read protien ids from proteinGroups.txt");
                            }
                            strProteinIdTemp = strLine.Substring(iBegain + 1, iEnd - iBegain - 1);  // Assume that the second column is "Majority protein IDs"
                            if(strProteinLastIdTemp!=strProteinIdTemp)
                            {
                                foreach (string item in ProteinsInFasta)
                                {
                                    if (item.IndexOf(strProteinIdTemp) != -1)
                                    {
                                        strProteinIdentifier = "";
                                        ProteinNames.Add(new ProteinName { ProteinHeader = item, ProteinIdentifier = strProteinIdentifier, ProteinIdentifierInInputFiles = strProteinIdTemp });
                                        strProteinLastIdTemp = strProteinIdTemp;
                                        iProteins++;
                                        break;
                                    }
                                }
                            }
                            
                            strLine = streamReader.ReadLine();
                        }
                        fasta_grid.DataContext = ProteinNames;
                        streamReader.Close();
                    }
                    catch (IOException ex)
                    {
                        System.Windows.MessageBox.Show("An IOException has been thrown!" + ex.Message + "\n Please Set the input directory in advance!");
                        return;
                    }
                    catch (ArgumentException ex)
                    {
                        System.Windows.MessageBox.Show("Cannot open " + ex.Message);
                        return;
                    }
                }
                else if(InputType=="SWATH")
                {
                    try
                    {
                        string strResultPath = getInputPathText.Invoke();
                        FileStream aFile = new FileStream(strResultPath, FileMode.Open);
                        StreamReader streamReader = new StreamReader(aFile);
                        ProteinNames = new ObservableCollection<ProteinName>();

                        iProteins = 0;
                        strLine = streamReader.ReadLine();
                        strLine = streamReader.ReadLine();
                        int iBegain, iEnd;
                        string strProteinIdTemp;
                        string strProteinLastIdTemp="";
                        while (strLine != null && iProteins < iShowLines)
                        {
                            while (strLine == "")
                            {
                                strLine = streamReader.ReadLine();
                            }

                            iBegain = 0;
                            iEnd = strLine.IndexOf(",", iBegain);
                            if (iEnd == -1)
                            {
                                System.Windows.MessageBox.Show("Cannot read protien ids from " + strResultPath);
                            }
                            strProteinIdTemp = strLine.Substring(iBegain, iEnd - iBegain);  // Assume that the first column is "Protein"
                            if(strProteinIdTemp!=strProteinLastIdTemp)
                            {
                                foreach (string item in ProteinsInFasta)
                                {
                                    if (item.IndexOf(strProteinIdTemp) != -1)
                                    {
                                        strProteinIdentifier = "";
                                        ProteinNames.Add(new ProteinName { ProteinHeader = item, ProteinIdentifier = strProteinIdentifier, ProteinIdentifierInInputFiles = strProteinIdTemp });
                                        strProteinLastIdTemp = strProteinIdTemp;
                                        iProteins++;
                                        break;
                                    }
                                }
                            }                            
                            strLine = streamReader.ReadLine();
                        }
                        fasta_grid.DataContext = ProteinNames;
                        streamReader.Close();
                    }
                    catch (IOException ex)
                    {
                        System.Windows.MessageBox.Show("An IOException has been thrown!" + ex.Message + "\n Please Check the input file path again.");
                        return;
                    }
                    catch (ArgumentException ex)
                    {
                        System.Windows.MessageBox.Show("Cannot open " + ex.Message);
                        return;
                    }
                }
                else if (InputType == "mzQuantML")
                { // Just load iShowLines proteins
                    try
                    {
                        ProteinNames = new ObservableCollection<ProteinName>();
                        int iShowProteins = 0;
                        string strProteinIdTemp;
                        if(iProteins<100)
                        {
                            iShowLines = iProteins;
                        }
                        else
                        {
                            iShowLines = 100;
                        }                                       
                        while (iShowProteins < iShowLines)
                        {
                            strProteinIdentifier = "";
                            strProteinIdTemp = "";
                            ProteinNames.Add(new ProteinName { ProteinHeader = ProteinsInFasta[iShowProteins], ProteinIdentifier = strProteinIdentifier, ProteinIdentifierInInputFiles = strProteinIdTemp });
                            iShowProteins++;
                        }
                        fasta_grid.DataContext = ProteinNames;
                    }
                    catch (IOException ex)
                    {
                        System.Windows.MessageBox.Show("An IOException has been thrown!" + ex.Message + "\n Please Set the input directory in advance!");
                        return;
                    }
                    catch (ArgumentException ex)
                    {
                        System.Windows.MessageBox.Show("Cannot open " + ex.Message);
                        return;
                    }                
                }
            }
        }
        private void Test_button_Click(object sender, RoutedEventArgs e)
        {
            string strLine;
            string strProteinIdentifier;
            List<string> ProteinsInFasta = new List<string>();
            bool bIfSetRECorrectly=true;
            Match match;
            int iProteins = 0;
            int iFastaProteinsNum = 999;
            try
            {
                if (!File.Exists(txtfasta.Text))
                {
                    System.Windows.MessageBox.Show("Cannot open the fasta file: " + txtfasta.Text);
                    return;
                }
                FileStream aFile = new FileStream(txtfasta.Text, FileMode.Open);
                StreamReader streamReader = new StreamReader(aFile);
                strLine = streamReader.ReadLine();
                iProteins = 0;
                while (strLine != null&&iProteins<=iFastaProteinsNum)
                {
                    while (strLine == "")
                    {
                        strLine = streamReader.ReadLine();
                    }
                    if (strLine[0] == '>')
                    {
                        ProteinsInFasta.Add(strLine);
                        iProteins++;
                        try
                        {
                            match = Regex.Match(strLine, txtRegularExpression.Text);
                            if (!match.Success)
                            {
                                System.Windows.MessageBox.Show("Cannot parse " + strLine + " by " + txtRegularExpression.Text);
                                return;
                            }
                        }
                        catch (ArgumentException)
                        {
                            System.Windows.MessageBox.Show("Seting a wrong regular expression.");
                            return;
                        }

                    }
                    strLine = streamReader.ReadLine();
                }
                fasta_status.Content ="Load protein headers from the fasta file.";
                streamReader.Close();
            }
            catch (ArgumentException){               
                System.Windows.MessageBox.Show("Please set the fasta file first!");
                return;
                }
            catch (IOException ex)
            {
                System.Windows.MessageBox.Show("An IOException has been thrown!" + ex.Message);
                return;
            }

            try
            {            
                string InputType = getInputType.Invoke();
                if (InputType == "maxquant")
                {
                    string strResultPath = getInputPathText.Invoke();
                    strResultPath = strResultPath + "\\proteinGroups.txt";
                    FileStream aFile = new FileStream(strResultPath, FileMode.Open);
                    StreamReader streamReader = new StreamReader(aFile);
                    ProteinNames = new ObservableCollection<ProteinName>();
                    int iShowLines = 100;
                    iProteins = 0;
                    strLine = streamReader.ReadLine();
                    strLine = streamReader.ReadLine();
                    int iBegain, iEnd;
                    string strProteinIdTemp;
                    string strProteinLastIdTemp = "";
                    bool bIfFindIndentifiedProteins = false;
                    while (strLine != null && iProteins < iShowLines)
                    {
                        while (strLine == "")
                        {
                            strLine = streamReader.ReadLine();
                        }
                        iBegain = strLine.IndexOf("\t");
                        if (iBegain == -1)
                        {
                            System.Windows.MessageBox.Show("Cannot read protien ids from proteinGroups.txt");
                            return;
                        }
                        iEnd = strLine.IndexOf("\t", iBegain + 1);
                        if (iEnd == -1)
                        {
                            System.Windows.MessageBox.Show("Cannot read protien ids from proteinGroups.txt");
                            return;
                        }
                        strProteinIdTemp = strLine.Substring(iBegain + 1, iEnd - iBegain - 1);  // Assume that the second column is "Majority protein IDs"
                       if(strProteinLastIdTemp!=strProteinIdTemp)
                       {
                           foreach (string item in ProteinsInFasta)
                           {
                               if (item.IndexOf(strProteinIdTemp) != -1)
                               {
                                   bIfFindIndentifiedProteins = true;
                                   match = Regex.Match(item, txtRegularExpression.Text);
                                   if (!match.Success)
                                   {
                                       System.Windows.MessageBox.Show("Cannot parse " + item + " by " + txtRegularExpression.Text);
                                       return;
                                   }
                                   else
                                   {
                                       strProteinIdentifier = match.Groups[1].Value;
                                       if (strProteinIdentifier != strProteinIdTemp)
                                       {
                                           bIfSetRECorrectly = false;
                                       }
                                       ProteinNames.Add(new ProteinName { ProteinHeader = item, ProteinIdentifier = strProteinIdentifier, ProteinIdentifierInInputFiles = strProteinIdTemp });
                                       strProteinLastIdTemp = strProteinIdTemp;
                                       iProteins++;
                                       break;
                                   }
                               }
                           } // for proteinsInFasta
                       }                       
                       strLine = streamReader.ReadLine();
                    }
                    fasta_grid.DataContext = ProteinNames;
                    streamReader.Close();
                    if(!bIfFindIndentifiedProteins)
                    {
                        System.Windows.MessageBox.Show("Cannot parse any identified protein from the fasta file! Please check the fasta file and the input files!");
                        return;
                    }
                    if (bIfSetRECorrectly)
                    {
                        System.Windows.MessageBox.Show("The regular expression is correct!");
                        return;
                    }
                    else
                    {
                        System.Windows.MessageBox.Show("Seting a wrong regular expression. \nPlease set it again to make sure " +
                            "that the identifiers extracted from the fasta file by the regular expression are the same as the protein IDs in input files.");
                    }
                }
                else if(InputType=="SWATH")
                {

                    string strResultPath = getInputPathText.Invoke();
                    FileStream aFile = new FileStream(strResultPath, FileMode.Open);
                    StreamReader streamReader = new StreamReader(aFile);
                    ProteinNames = new ObservableCollection<ProteinName>();
                    int iShowLines = 100;
                    iProteins = 0;
                    strLine = streamReader.ReadLine();
                    strLine = streamReader.ReadLine();
                    int iBegain, iEnd;
                    string strProteinIdTemp;
                    string strProteinLastIdTemp = "";
                    bool bIfFindIndentifiedProteins = false;
                    while (strLine != null&&iProteins < iShowLines)
                    {
                        while (strLine == "")
                        {
                            strLine = streamReader.ReadLine();
                        }

                        iBegain = 0;
                        iEnd = strLine.IndexOf(",", iBegain);
                        if (iEnd == -1)
                        {
                            System.Windows.MessageBox.Show("Cannot read protien ids from " + strResultPath);
                            return;
                        }
                        strProteinIdTemp = strLine.Substring(iBegain, iEnd - iBegain);  // Assume that the first column is "Protein"
                        if(strProteinLastIdTemp!=strProteinIdTemp)
                        {
                            foreach (string item in ProteinsInFasta)
                            {
                                if (item.IndexOf(strProteinIdTemp) != -1)
                                {
                                    bIfFindIndentifiedProteins = true;
                                    match = Regex.Match(item, txtRegularExpression.Text);
                                    if (!match.Success)
                                    {
                                        System.Windows.MessageBox.Show("Cannot parse " + item + " by " + txtRegularExpression.Text);
                                        return;
                                    }
                                    else
                                    {
                                        strProteinIdentifier = match.Groups[1].Value;
                                        if (strProteinIdentifier != strProteinIdTemp)
                                        {
                                            bIfSetRECorrectly = false;
                                        }
                                        ProteinNames.Add(new ProteinName { ProteinHeader = item, ProteinIdentifier = strProteinIdentifier, ProteinIdentifierInInputFiles = strProteinIdTemp });
                                        strProteinLastIdTemp = strProteinIdTemp;
                                        iProteins++;
                                        break;
                                    }
                                }
                            } // for proteinsInFasta
                        }
                        
                        strLine = streamReader.ReadLine();
                    }
                    fasta_grid.DataContext = ProteinNames;
                    streamReader.Close();
                    if(!bIfFindIndentifiedProteins)
                    {
                        System.Windows.MessageBox.Show("Cannot parse any identified protein from the fasta file! Please check the fasta file and the input files!");
                        return;
                    }
                    if (bIfSetRECorrectly)
                    {
                        System.Windows.MessageBox.Show("The regular expression is correct!");
                        return;
                    }
                    else
                    {
                        System.Windows.MessageBox.Show("Seting a wrong regular expression. \nPlease set it again to make sure " +
                            "that the identifiers extracted from the fasta file by the regular expression are the same as the protein IDs in input files.");
                    }
                }
                else if(InputType=="mzQuantML")
                { // Just show iShowLines proteins

                    int iShowLines;
                    if(iProteins<100)
                    {
                        iShowLines = iProteins;
                    }
                    else
                    {
                        iShowLines = 100;
                    }

                    int iShowProteins = 0;
                    string strProteinIdTemp;
                    string strLastProteinIdentifier="";
                    ProteinNames = new ObservableCollection<ProteinName>();
                    while (iShowProteins < iShowLines)
                    {
                        match = Regex.Match(ProteinsInFasta[iShowProteins], txtRegularExpression.Text);
                        if (!match.Success)
                        {
                            System.Windows.MessageBox.Show("Cannot parse " + ProteinsInFasta[iShowProteins] + " by " + txtRegularExpression.Text);
                            return;
                        }
                        else
                        {
                            strProteinIdentifier = match.Groups[1].Value;
                            if(strLastProteinIdentifier!=strProteinIdentifier)
                            {
                                strProteinIdTemp = "";
                                ProteinNames.Add(new ProteinName { ProteinHeader = ProteinsInFasta[iShowProteins], ProteinIdentifier = strProteinIdentifier, ProteinIdentifierInInputFiles = strProteinIdTemp });
                                strLastProteinIdentifier = strProteinIdentifier;
                                iShowProteins++;
                            }
                        }
                    }
                        
                    }
                    fasta_grid.DataContext = ProteinNames;                   
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
        }

        private void Parse_OK_button_Click(object sender, RoutedEventArgs e)
        {
            if(ChangeTextEvent!=null)
            {
                ChangeTextEvent.Invoke(txtRegularExpression.Text);
                this.Close();
            }
        }

        private void txtRegularExpression_DataContextChanged(object sender, DependencyPropertyChangedEventArgs e)
        {

        }

        private void txtRegularExpression_TextChanged(object sender, TextChangedEventArgs e)
        {
        }
    }
   public class RegularExpressionExample:INotifyPropertyChanged
    {

       public event PropertyChangedEventHandler PropertyChanged;
       private string m_header;
       private string m_re;
       private string m_Identifier;
       public string header
       {
           get { return m_header; }
           set
           {
               m_header = value;
               if (PropertyChanged != null)
               {
                   PropertyChanged(this, new PropertyChangedEventArgs("header"));                
               }
           }
       }
       public string re
       {
           get { return m_re; }
           set
           {
               m_re = value;
               if (PropertyChanged != null)
               {
                   PropertyChanged(this, new PropertyChangedEventArgs("re"));
               }
           }
       }
       public string Identifier
       {
           get { return m_Identifier; }
           set
           {
               m_Identifier = value;
               if (PropertyChanged != null)
               {
                   PropertyChanged(this, new PropertyChangedEventArgs("Identifier"));
               }
           }
       }
    }

    public class ProteinName:INotifyPropertyChanged
    {
        private string m_ProteinHeader;
        private string m_ProteinIdentifier;
        private string m_ProteinIdentifierInInputFiles;
        public event PropertyChangedEventHandler PropertyChanged;
        public string ProteinHeader 
        {
            get
            {
                return m_ProteinHeader;
            }
            set
            {
                m_ProteinHeader = value;
                if(PropertyChanged!=null)
                {
                    PropertyChanged(this, new PropertyChangedEventArgs("ProteinHeader"));
                }
            }
        }
        public string ProteinIdentifier
        {
            get
            {
                return m_ProteinIdentifier;
            }
            set
            {
                m_ProteinIdentifier = value;
                if (PropertyChanged != null)
                {
                    PropertyChanged(this, new PropertyChangedEventArgs("ProteinIdentifier"));
                }
            }
        }
        public string ProteinIdentifierInInputFiles
        {
            get
            {
                return m_ProteinIdentifierInInputFiles;
            }
            set
            {
                m_ProteinIdentifierInInputFiles = value;
                if (PropertyChanged != null)
                {
                    PropertyChanged(this, new PropertyChangedEventArgs("ProteinIdentifierInInputFiles"));
                }
            }
        }
    }


}
