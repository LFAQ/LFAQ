   M S B A M L   `   `   `  ����   LFAQ    LFAQ.MainWindow   \ XPresentationFramework, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35R NWindowsBase, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35W SPresentationCore, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35D 9http://schemas.microsoft.com/winfx/2006/xaml/presentation    R NSystem.Xaml, Version=4.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089:x,http://schemas.microsoft.com/winfx/2006/xaml     -     �Name$  window��$	��500��$	��700�� � Title$
 LFAQ�� � Icon$/ )pack://application:,,,/Resources/logo.png��.��]� -    ]�Name$
 dock��.+�o� $	��Top=�
�������
���   �$
��Left=�.@�m� -   $��MenuFile��$��_File��.@�m� -   $��Load parameters��m� -   $��Save parameters���� m� -   $��_Exit�� � ToolTip$ Exit the application��m� -   $��menuHelp��$��_Help��.@�m� -   $��About��m� -	   $��Online User Guide���� -
   $��	statusBar��$��20��$��Bottom=�.@��� -    ��Name$ statBarText��$
��Wrap=�$��Ready�� J� OpacityMask �	#FFCA3131� -    �Name$ grid1��$	��Top=�.+��� -    ��Name$ 
Global_tab��$��0,0,-0.4,-0.6q����� *ActualWidth+	 �� ElementName$	 grid1��
 �� Mode$
 OneWay=����� *ActualHeight+$	 grid1��$
 OneWay=�.@��� $
��Data��.��� ���	#FFE5E5E5$��	0,2,0.2,0q� � RowDefinitions �� $��* ��� $��* ��� $��* ��� $��3* ��� $��* ��� $��* �.+��� $��0�����	#FFFCF8F8$��0.5,0.5,0.5,0.5q�._��� -    ��Name$ stackPanel0��$��0�� �� Orientation$ 
Horizontal=�$��
0,0,-0.2,0q�.+��� $��
Input type��$��Center=� ^� ShowDuration$
 2000�� �� .+��� ;9Please choose the software used to identify the proteins.�� -   $��MaxQuantCheckBox��$��MaxQuant��$
��Left=�$��Center=� d� 	IsChecked$
 TrueY��� -   $��SWATHCheckBox��$��	SWATH 2.0��$
��Left=�$��20,0,0,0q�$��Center=��� -   $��mzQuantMLCheckBox��$��	mzQuantML��$
��Left=�$��20,0,0,0q�$��Center=��� $��1�����	#FFFCF8F8$��0.5,0.5,0.5,0.5q�._��� -   $��stackPanel1��$��1��$ 
Horizontal=�$��
0,0,-0.2,0q�.+��� $��
Fasta file��$
��Left=�$��Center=�$
 2000�� �� .+��� .,The file path of proteins database (*.fasta)�� -    ��Name$ txtfasta��$��72,0,0,0q�$��Center=�$��Center=�$	��150��$	��400�� �� TextWrapping$
 Wrap=��� -    ��Name$ fasta_browse��$��Browse��$��10,0,0,0q�$��Center=�$��Center=��� $��2�����	#FFFCF8F8$��0.5,0.5,0.5,0.5q�._��� -   $��stackPanel2��$��2��$ 
Horizontal=�$��
0,0,-0.2,0q�.+��� $��Identifier parsing rule��$
��Left=�$��Center=� �� .+��� $
��Wrap=�QOThe regular expression used to extract protein identifiers from the fasta file.�� -   $ txtfastaType��$��Center=�$��23��$
 Wrap=�$��Center=�$	��150��$��>(.*?)\s���� -   $ 
fasta_type��$	��Set��$��10,0,0,0q�$��Center=�$��Center=����� *ActualWidth+$	 fasta_browse��$
 OneWay=��� $��3�����	#FFFCF8F8$��0.5,0.5,0.5,0.5q�._�� -   $��
GridPanel3�� � ColumnDefinitions �� $��* ��� $��* �.+��� $��0��$��Center=�.+��� $9��3Does the input protein list contain decoy proteins?���� $ 
Horizontal=�$��5,0,0,0q�.+��� -   $��ContainDecoyProtein��$	��Yes���� -   $��NoDecoyProtein��$��No��$��20,0,0,0q�$
 TrueY��� $��0,15,0,0q�.+��� $9��3If yes, please provide the prefix of decoy proteins���� -   $��txtDecoyPrefix��$��5,0,0,0q�$	��100��$
��Left=�.���� $��1��$��Center=�.+��� $?��9Does the input protein list contain contaminant proteins?���� $ 
Horizontal=�$��5,0,0,0q�.+��� -   $��ExistContaProteins��$	��Yes���� -   $��NoContaProteins��$��No��$��20,0,0,0q�$
 TrueY��� $��0,15,0,0q�.+��� $?��9If yes, please provide the prefix of contaminant proteins���� -   $��txtContaminantPrefix��$��5,0,0,0q�$	��100��$
��Left=�.���� $��4�����	#FFFCF8F8$��0.5,0.5,0.5,0.5q�._�� .+��� -   $��InputDirectoryStack��$ 
Horizontal=�.+��� $��Input directory��$��Center=��� -    $ txtInputDirectoryPath��$��	30,0,10,0q�$	��400��$	��550��$
 Wrap=�$��Center=�$��Center=��� -!   $ InputDirectory_Browse��$��Browse��$
��Left=�$��Center=�$	��600���� -"   $��InputFileStack��$ 
Horizontal=�$y�Hidden=�.+��� $��
Input file��$��Center=��� -#   $ txtInputFilePath��$��	30,0,10,0q�$	��300��$	��350��$
 Wrap=�$��Center=�$��Center=��� -$   $ InputFile_Browse��$��Browse��$
��Left=�$��Center=�$	��600��� $��5��.+��� -%   $��Run_buttonInData��$��1��$��Run LFAQ��$��Center=�$��Center=�$��42��$��73���� $��Quantification��.��� ���	#FFE5E5E5 �� $��8* ��� $��4* ��� $��2* �.+��� -&   $��RegressionParamGroupBox��$��0��$��
Regression��.���  �� $��* ��� $��* �.+�� $��0��.+��� $��Regression method��$
��Left=�$��10,0,0,0q�$	��Top=�$	��125�� �� .+��� <:Please choose the regression method for Q-factor learning.�� -'   $��ComboBoxRegressionMethod��$
��Left=�$��	140,0,0,0q�$	��Top=�$	��120�� �� SelectedIndex$ 1��.@��� -(   $��CbiStepwise��$��stepwise��$
��Left=��� -)   $��CbiBART��$
��BART��$
��Left=�� -*   $ stepwisePanel��$��0��$��0,10,0,0q�$	��474��$y�Hidden=� �� $
��auto ��� $��* ��� $
��auto ��� $��* � �� $��* ��� $��* �+��� $��0��$��0��$��alpha1��$��0,5,2,0q�$��Center=�$��Bottom=� �� .+��� :8The alpha1 should be a numerical number between 0 and 1.�� $��1��$��0��$��alpha2��$��0,5,2,0q�$��Center=�$	��Top=� �� .+��� :8The alpha2 should be a numerical number between 0 and 1.�� -+   $ 	txtAlpha1��$��0��$��1��$��68��$��32,5,0,5q�$
 Wrap=�$
��0.95��$��Bottom=��� -,   $��	txtAlpha2��$��1��$��1��$
 Wrap=�$��32,5,0,0q�$��68��$
��0.95��$	��Top=�� --   $ 	BARTPanel��$��0��$��0,30,0,0q�$y�Visible=� �� $��* ��� $��* �+�� $��0�� �� $��* ��� $��* �.+��� $��0��$ 
Horizontal=�.+��� $��Alpha��$��Center=�$��Center=� �� .+��� =;The base parameter for the tree prior, ranging from 0 to 1.�� -.   $ txtAlpha��$
 Wrap=�$��	144,0,0,0q�$
��0.85��$��68��$
��Auto��$��Center=�.���� $��1��$ 
Horizontal=�.+��� $
��Beta��$��Center=�$	��Top=� �� .+��� NLThe power parameter for the tree prior, ranging from 0 to positive infinite.�� -/   $ txtBeta��$
 Wrap=�$��	151,0,0,0q�$��68��$	��1.6��$	��Top=�� $��1�� �� $��* ��� $��* �.+��� $��0��$ 
Horizontal=�.+��� $��k��$��3q�$��Center=�$��Center=� �� .+��� SQThe number of standard deviations of the dependent variables in the training set.�� -0   $
 txtK��$
 Wrap=�$��	236,0,0,0q�$��68��$��2��$��Center=��� $��1��$ 
Horizontal=�.+��� $��Number of trees��$��Center=�$	��Top=� �� .+��� +)The number of trees to train in the BART.�� -1   $ txtNumberOfTrees��$
 Wrap=�$	��200��$��	154,0,0,0q�$��68��$	��Top=�� $��1�� �� $��* ��� $��* �.+�� $��0�� �� $��* ��� $��* ��� $��* ��� $��* �.+��� $��0��$ 
Horizontal=�.+��� $��Max missed cleavage(s)��$
��Left=�$��Center=� �� .+��� SQThe maximum number of missed cleavages of a peptide in the theoretical digestion.�� -2   $ txtMaxMissedCleave��$��42,0,0,0q�$
 Wrap=�$��Center=�$��68��$��0���� $��1��$ 
Horizontal=�.+��� $!��Shortest length of peptides��$��Center=� �� .+��� HFThe allowed shortest length of a peptide in the theoretical digestion.�� -3   $ txtPepShotest��$��20,0,0,0q�$
 Wrap=�$��Center=�$��68��$��6���� $��2��$ 
Horizontal=�.+��� $ ��Longest length of peptides��$��Center=� �� .+��� GEThe allowed longest length of a peptide in the theoretical digestion.�� -4   $ txtPepLongest��$��Right=�$��22,0,0,0q�$
 Wrap=�$��Center=�$��68��$��30���� $��3��$ 
Horizontal=�.+��� $��Enzyme��$
��Left=�$��Center=� �� .+��� ,*The enzyme used for theoretical digestion.�� -5   $ 	txtEnzyme��$��Right=�$��	132,0,0,0q�$
 Wrap=�$��Center=�$��68��$��trypsin��� $��1�� �� $��* ��� $��2* ��� $��2* ��� $��2* �.+��� $��1��$ 
Horizontal=�.+��� $��0,2,0,0q�$/��)Does the sample contain standard proteins���� -6   $��ComboBoxContainStand��$��22,0,0,0q�$ 0��$��20��.@��� -7   $��CbiYes��$	��yes��$
��Left=��� -8   $��CbiNo��$��no��$
��Left=�� -9   $��StandProteinGrid��$��2��$��2��$y�Visible=� �� $��* ��� $��* �.+�� $��0��.+��� $@��:If yes, please provide the identifier of standard proteins��� $��1��.+��� -:   $��txtIdentifierOfStand��$��5,0,0,0q�$	��200��$	��100��$
��Left=�$��Center=�$	��ups��.���� -;   $ SetStandProteins��$��Set standard proteins��$��	195,0,0,0q�$
��Left=�$��Center=��� -<   $��QuantParamGroupBox��$��1��$��Quantification��.��� .+��� -=   $ 
stackPanel��$ 
Horizontal=�.+��  �� $��* ��� $��* �.+��� ->   $��IfCalculateiBAQ��$��0��$��Calculate iBAQ��$��Bottom=�$��0,0,50,0q�$
 TrueY� �� .+��� trIf MaxQuant result is used as input, iBAQ values come from MaxQuant result. Otherwise, iBAQ is calculated by LFAQ.�� -?   $��IfCalculateTop3��$��1��$��Calculate Top3��$��Center=�$��0,0,50,0q�$
 TrueY��� $��Result directory��$��Center=�$��0,20,0,0q� �� .+��� MKPlease choose the directory where the quantification results will be saved.�� -@   $ txtResultPath��$	��300��$	��360��$
 Wrap=�$��Center=�$��0,20,0,0q��� -A   $ Result_Browse��$��Browse��$��Center=�$��	10,20,0,0q�� $��2��.+��� -B   $ ��Run_ButtonInQuantification��$��1��$��Run LFAQ��$��Center=�$��Center=�$��42��$��73���� -C    ��Name$ VisualizationTabItem��$��Visualization��$��20��$	��Top=�$��-1.8,-2,1.6,0q�.��� -D   $��VisualizationGrid�� �� $
��618* ��� $
��382* �.+�� $��0����.  �� $
��618* ��� $
��382* �.+�&   System.Windows.Controls.DataGrid  -E   $��LFAQ_result_grid��$��0��   AutoGenerateColumns .  $   HorizontalScrollBarVisibility$ Visible=�$ 0������    Columns 0  *System.Windows.Controls.DataGridTextColumn  ,  &System.Windows.Controls.DataGridColumn   Header$ Major Protein IDs��   Width5  /System.Windows.Controls.DataGridLengthConverter$ *    MaxWidth$	 200Q�   MinWidth$ 50Q�   FontSize$ 8Q�1  +System.Windows.Controls.DataGridBoundColumn   Binding �� *strProteinName+   HeaderStyle ��   �� 
TargetType=  7System.Windows.Controls.Primitives.DataGridColumnHeader"   ��� ! �� Property! � ��" �� Value$" 8(�  $ Protein IDs ��$ * $	 200Q�$ 50Q�$ 8Q� �� *strProteinIDs+ �� "   ��� ! � ��$" 8(�  $ Experiment Name��$ 50Q�$ 8Q� �� *strExperiment+ �� "   ��� ! � ��$" 8(�  $
 LFAQ��$ 30Q�$ 8Q� �� *	LFAQpep+ �� "   ��� ! � ��$" 8(�  $
 iBAQ��$ 30Q�$ 8Q� �� *MaxquantIBAQ+ �� "   ��� ! � ��$" 8(�  $
 Top3��$ 30Q�$ 8Q� �� *Top3+ �� "   ��� ! � ��$" 8(�  $ PredictedMol(LFAQ)��$ * $ 60Q�$ 8Q� �� *PredictedMolOfLFAQpep+ �� "   ��� ! � ��$" 8(�  $ PredictedMol(iBAQ)��$ * $ 60Q�$ 8Q� �� *PredictedMolOfiBAQ+ �� "   ��� ! � ��$" 8(�  $ PredictedMol(Top3)��$ * $ 60Q�$ 8Q� �� *PredictedMolOfTop3+ �� "   ��� ! � ��$" 8(��� -F   $��canvasOfOneProtein��$��1��$��5,20,5,5q����	#FFFBF6F6.+��� -G   $��OneProtein_image��$
��Fill=����� *ActualWidth+$	 canvasOfOneProtein��$
 OneWay=����� *ActualHeight+$	 canvasOfOneProtein��$
 OneWay=�� $��1�� �� $
��618* ��� $
��382* �.+�  -H   $��dataGrid_StandProteins��$��1�� .  ����    $ Major Protein IDs��$ * $	 200Q�$ 50Q�$ 8Q� �� *strProteinName+ �� "   ��� ! � ��$" 8(�  $ Protein IDs��$ * $	 200Q�$ 50Q�$ 8Q� �� *strProteinIDs+ �� "   ��� ! � ��$" 8(�  $ Experiment Name��$ 50Q�$ 8Q� �� *strExperiment+ �� "   ��� ! � ��$" 8(�  $ SpikedInMols��$ * $ 30Q�$ 8Q� �� *SpikedInMols+ �� "   ��� ! � ��$" 8(�  $
 LFAQ��$ 30Q�$ 8Q� �� *	LFAQpep+ �� "   ��� ! � ��$" 8(�  $
 iBAQ��$ 30Q�$ 8Q� �� *MaxquantIBAQ+ �� "   ��� ! � ��$" 8(�  $
 Top3��$ 30Q�$ 8Q� �� *Top3+ �� "   ��� ! � ��$" 8(�  $ PredictedMol(LFAQ)��$ * $ 60Q�$ 8Q� �� *PredictedMolOfLFAQpep+ �� "   ��� ! � ��$" 8(�  $ PredictedMol(iBAQ)��$ * $ 60Q�$ 8Q� �� *PredictedMolOfiBAQ+ �� "   ��� ! � ��$" 8(�  $ PredictedMol(Top3)��$ * $ 60Q�$ 8Q� �� *PredictedMolOfTop3+ �� "   ��� ! � ��$" 8(�� $��1�� �� $��2* ��� $��8* �.+�  -I   $��dataGrid_SelectExperiment��$��1�� .  $��	5,20,0,10q�#   
IsReadOnly# . ����    $ 
Experiment��$ * $ 50Q�$ 8Q� �� *strExperimentName+ �� "   ��� ! � ��$" 8(��� -J   $ ��Name$$ canvas��$��1��$��5,20,5,5q�.+��� -K   $��
LFAQ_image��$
��Fill=����� *ActualWidth+$	 canvas��$
 OneWay=����� *ActualHeight+$	 canvas��$
 OneWay=�