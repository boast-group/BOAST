function scanner_param = SetDefaultScannerParam

our_scanner.name = 'Trio';
our_scanner.B0 = 3;
% well this is not directly a scanner parameter, however it depends on the
% field strength and therefore scanner
our_scanner.T2s = 45*10^-3; 	% T2* at 3 T (Wansapura et al., JMRI 1999)

scanner_param = our_scanner;

end
