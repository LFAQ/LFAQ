   M S B A M L   `   `   `  ����	  SPEAQ     SPEAQ.SetRegularExpression   5   	   \ XPresentationFramework, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35R NWindowsBase, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35W SPresentationCore, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35D 9http://schemas.microsoft.com/winfx/2006/xaml/presentation    5   	   R NSystem.Xaml, Version=4.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089:x,http://schemas.microsoft.com/winfx/2006/xaml     5   	     � Title$  Set regular expression��5   	   .��5      �  � RowDefinitions 5   
   �� 5      $��* �6   6*   �� 5      $��* �6   6)   �� 5	      $��6* �6   6*   �� 5
      $��2* �6   6*   .+�5   
   ]� $��0��6   .+�5      �� $��Identifier parse rule��6   $
��Left=�65   $��Center=�6Q   �� 5      -   5        $��txtRegularExpression��5      $��30,0,0,0q�62   $
��Left=�6D   $��23��6_    �� TextWrapping$
 Wrap=�6l   $��Center=�6�   $	��120��6�   6�   5      ]� 5   
   $��1��6   .+�5      �� $��
Fasta file��6   $
��Left=�6*   $��Center=�6F   �� 5      -   5         ��Name$ txtfasta��5      $��30,0,0,0q�6(   $
��Left=�6:   $
 Wrap=�6V   $��Center=�6k   $��80��6�   $	��280��6�   �� 5      -   5         ��Name$ fasta_browse��5      $��Browse��6+   $��10,0,0,0q�6<   $��Center=�6N   $��Center=�6k   6�   �� 5      -   5        $ Test_button��5       ^� ToolTip$� �Note that the identifiers extracted from the fasta file by the regular expression should be the same as the protein IDs in input files.��6*    ^� InitialShowDelay$ 5��5       ^� 	Placement$ Bottom=�69   $
��Test��6[   $��
0,3,10,3.2q�6j   $
��Left=�6~   $��Center=�6�   6�   �� 5      -   5        $ Parse_OK_button��5      $��OK��6.   $��1,0,10,0q�6;   $
��Left=�6M   $��Center=�6h   6�   5      ]� 5   
   $��2��6   ��. .+�5      �� -   5         ��Name$ fasta_status��5      $	��Top=�6*   5      &   System.Windows.Controls.DataGrid  5      -   5        $��
fasta_grid��5      	   AutoGenerateColumns	 .  ��6E   �� 
   Columns
 5      0  *System.Windows.Controls.DataGridTextColumn  5      ,  &System.Windows.Controls.DataGridColumn   Header$ Protein Header��6L      MaxWidth$	 500Q�6d   1  +System.Windows.Controls.DataGridBoundColumn   Binding 6)   �� *ProteinHeader+6d     5      $ 
Identifier��6P      Width5  /System.Windows.Controls.DataGridLengthConverter$ * 6d    6)   �� *ProteinIdentifier+6d     5       $ Identifier In Input Files��6\   $ * 6    6)   �� *ProteinIdentifierInInputFiles+6   5"      5#      ]� 5$   
   ��. $��3��6)   .+�5%      �� $��	Examples:��6   5&        5'      -   5        $��examples��5'      	 .  ��6'   �� 
 5(        5)      $ Protein Header��6D   $	 300Q�6\    6)   �� *header+6\     5*      $ Regular Expression��6@   $ * 6\    6)   �� *re+6\     5+      $ 
Identifier��6H   $ * 6\    6)   �� *
Identifier+6\   5-      5.      50      51      