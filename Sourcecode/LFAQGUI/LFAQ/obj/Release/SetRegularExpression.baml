   M S B A M L   `   `   `  ����   LFAQ    LFAQ.SetRegularExpression   \ XPresentationFramework, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35R NWindowsBase, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35W SPresentationCore, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35D 9http://schemas.microsoft.com/winfx/2006/xaml/presentation    R NSystem.Xaml, Version=4.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089:x,http://schemas.microsoft.com/winfx/2006/xaml       � Title$  Set regular expression��.���  � RowDefinitions �� $��* ��� $��* ��� $��6* ��� $��2* �.+�]� $��0��.+��� $��
Fasta file��$
��Left=�$��Center=��� -    ��Name$ txtfasta��$��92,0,0,0q�$
��Left=� �� TextWrapping$
 Wrap=�$��Center=�$	��300��$	��500���� -    ��Name$ fasta_browse��$��Browse��$��10,0,0,0q�$
��Left=�$��Center=�]� $��1��.+��� $��Identifier parsing rule��$
��Left=�$��Center=��� -   $��txtRegularExpression��$��20,0,0,0q�$
��Left=�$��23��$
 Wrap=�$��Center=�$	��120��$��>(.*?)\s���� -   $ Test_button�� ^� ToolTip$� �Note that the identifiers extracted from the fasta file by the regular expression should be the same as the protein IDs in the input files�� ^� InitialShowDelay$ 5�� ^� 	Placement$ Bottom=�$
��Test��$��10,0,0,0q�$
��Left=�$��Center=��� -   $ Parse_OK_button��$��OK��$��	155,0,0,0q�$
��Left=�$��Center=�]� $��2����. .+��� -    ��Name$ fasta_status��$	��Top=�&   System.Windows.Controls.DataGrid  -   $��
fasta_grid��	   AutoGenerateColumns	 .  ���� 
   Columns
 0  *System.Windows.Controls.DataGridTextColumn  ,  &System.Windows.Controls.DataGridColumn   Header$ Protein header��   MaxWidth$	 500Q�1  +System.Windows.Controls.DataGridBoundColumn   Binding �� *ProteinHeader+  $ 
Identifier��   Width5  /System.Windows.Controls.DataGridLengthConverter$ *  �� *ProteinIdentifier+  $ Identifier in input files��$ *  �� *ProteinIdentifierInInputFiles+]� ��. $��3��.+��� $��	Examples:��  -   $��examples��	 .  ���� 
   $ Protein header��$	 300Q� �� *header+  $ Regular expression��$ *  �� *re+  $ 
Identifier��$ *  �� *
Identifier+