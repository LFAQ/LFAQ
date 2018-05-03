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

namespace LFAQ
{
    /// <summary>
    /// Interaction logic for SetStandProteins.xaml 
    /// </summary>
    //public  delegate string GetTextBoxContentHandler();
    public delegate string GetInputPathHandler();
    public delegate string GetResultPathHandler();
   // public delegate string GetInputTypeHandle();
    public partial class SetStandProteins : Window
    {
        public GetTextBoxContentHandler GetIdentifierOfStand;
        public GetInputPathHandler GetInputPath;
        public GetInputTypeHandle GetInputType;
        public GetResultPathHandler GetResultPath;
        public ObservableCollection<StandProteinAmount> StandProteinsAmounts { get; set; }
        public SetStandProteins()
        {
            InitializeComponent();
            try
            {

                StandProteinsAmounts = new ObservableCollection<StandProteinAmount>();
                StandProteinsAmounts.Clear();
                string strLine;
                int iStart, iEnd;
                string CurrentDirectory = System.Environment.CurrentDirectory;
                string StandProteinsPath = "StandardProteins.txt";
                string strProteinName;
                double dSpikedInOfProtein;
                using (StreamReader streamReader = new StreamReader(StandProteinsPath, Encoding.Default))
                {
                    strLine = streamReader.ReadLine();
                    strLine = streamReader.ReadLine();
                    while (strLine != null)
                    {
                        iStart = 0;
                        iEnd = strLine.IndexOf("\t");
                        strProteinName = strLine.Substring(iStart, iEnd - iStart);  //uniProt Accession Number

                        iStart = iEnd + 1;
                        iEnd = strLine.IndexOf("\t", iStart);   //Amount (fmol)
                        if (iEnd == -1)
                        {
                            iEnd = strLine.IndexOf("\n", iStart);
                        }
                        if (iEnd == -1)
                        {
                            dSpikedInOfProtein = Convert.ToDouble(strLine.Substring(iStart, strLine.Length - iStart));
                        }
                        else
                        {
                            dSpikedInOfProtein = Convert.ToDouble(strLine.Substring(iStart, iEnd - iStart));
                        }
                        

                        try
                        {
                            StandProteinsAmounts.Add(new StandProteinAmount { strProteinName = strProteinName, dAmount = dSpikedInOfProtein });
                        }
                        catch (System.Exception excep)
                        {
                            System.Windows.MessageBox.Show(excep.ToString());
                        }

                        strLine = streamReader.ReadLine();
                    }
                }
                StandProteins_grid.DataContext = StandProteinsAmounts;
            }
            catch (IOException ex)
            {
                System.Windows.MessageBox.Show("An IOException has been thrown!"+ex.Message);
                return;
            }
        }

        private void SaveStandProteins_click(object sender, RoutedEventArgs e)
        {
             string strIdentifierOfStand = GetIdentifierOfStand.Invoke();
            // Make sure all standard protein names contain the identifier;
            // and the the number of standard proetins with different amount should be greater than 3;
            Dictionary<double, string> dictStandProAmountAndNames=new Dictionary<double,string>();
            foreach (var item in StandProteinsAmounts)
            {
           
                if(item.strProteinName.IndexOf(strIdentifierOfStand)==-1)
                {
                    System.Windows.MessageBox.Show(strIdentifierOfStand + " is not the identifier of " + item.strProteinName);
                    return;
                }
            }
            string strIdentificationResultPath = GetInputPath.Invoke();
            string strQuantResultPath = GetResultPath.Invoke();
            string StandProteinsPath = strQuantResultPath + "\\StandardProteins.txt";

            // make sure that there are standard proteins (users setted) exist in the identification result.
            try
            {
                bool bIfStandAppeared = false;
                string strIdentifiedProteinsPath;
                string strInputType = GetInputType.Invoke();
                List<string> listIdentifiedProteinNames = new List<string>();
                if (strInputType == "maxquant")
                {
                    strIdentifiedProteinsPath = strIdentificationResultPath + "\\proteinGroups.txt";

                    FileStream aFile = new FileStream(strIdentifiedProteinsPath, FileMode.Open);
                    StreamReader streamReader = new StreamReader(aFile);
                    ObservableCollection<ProteinName> ProteinNames = new ObservableCollection<ProteinName>();
                    string strLine;
                    int iBegain, iEnd;
                    string strProteinIdTemp;

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
                            System.Windows.MessageBox.Show("Cannot read protien ids from proteinGroups.txt");
                        }
                        iEnd = strLine.IndexOf("\t", iBegain + 1);
                        if (iEnd == -1)
                        {
                            System.Windows.MessageBox.Show("Cannot read protien ids from proteinGroups.txt");
                        }
                        strProteinIdTemp = strLine.Substring(iBegain + 1, iEnd - iBegain - 1);  // Assume that the second column is "Majority protein IDs"
                        listIdentifiedProteinNames.Add(strProteinIdTemp);

                        foreach (var item in StandProteinsAmounts)
                        {
                           if (strProteinIdTemp.IndexOf(item.strProteinName) != -1)
                           // if (strProteinIdTemp==item.strProteinName)
                            {
                                bIfStandAppeared = true;
                                if (dictStandProAmountAndNames.ContainsKey(item.dAmount) == false)
                                {
                                    dictStandProAmountAndNames.Add(item.dAmount, item.strProteinName);
                                }
                                break;
                            }
                        }


                        strLine = streamReader.ReadLine();
                    }

                    streamReader.Close();
                    if (bIfStandAppeared == false)
                    {
                        System.Windows.MessageBox.Show("Cannot find any standard protein in input files.");
                        return;
                    }
                    
                    if (dictStandProAmountAndNames.Count < 3)
                    {
                        System.Windows.MessageBox.Show("There are not enough standard proteins identified in the input file. At least three standard proteins with different spiked-in amounts are needed!");
                        return;
                    }
                }
                else if (strInputType=="mzQuantML")
                {
                    //ToDo: Test if there are at least three standard proteins with different spiked-in amounts in input file
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

            SaveStandProteinsAmount(StandProteinsAmounts, StandProteinsPath);
            this.Close();
        }

        private void SaveStandProteinsAmount(ObservableCollection<StandProteinAmount> StandProteinsAmounts,string path)
        {
            try
            {
                StreamWriter writer = File.CreateText(path);
                writer.WriteLine("Protein Id\tAmount\t");
                foreach (var item in StandProteinsAmounts)
                {
                   writer.WriteLine(item.strProteinName + "\t" + item.dAmount+"\t");
                }
                writer.Close();
            }
            catch (IOException ex)
            {
                System.Windows.MessageBox.Show("An IOException has been thrown!"+ex.Message);
                return;
            }
        }
        public class StandProteinAmount
        {
            public event PropertyChangedEventHandler PropertyChanged;
            private string m_strProteinName;
            private double m_dAmount;

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
            public double dAmount
            {
                get { return m_dAmount; }
                set
                {
                    m_dAmount = value;
                    if (PropertyChanged != null)
                    {
                        PropertyChanged(this, new PropertyChangedEventArgs("dAmount"));
                    }
                }
            }
        }

    }
}
