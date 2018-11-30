try
    load SourceData;
catch
    SourceData = containers.Map;
end

SourceData = [ SourceData; ObtainSingleBEASpreadsheet( 'https://apps.bea.gov/national/FixedAssets/Release/XLS/Section2All_xls.xlsx', 'FAAt', '-A', { '201', '202' } ) ];
SourceData = [ SourceData; ObtainSingleBEASpreadsheet( 'https://apps.bea.gov/national/Release/XLS/Survey/Section1All_xls.xlsx', 'T', '-A', { '10103', '10105' } ) ];
SourceData = [ SourceData; ObtainSingleBEASpreadsheet( 'https://apps.bea.gov/national/Release/XLS/Survey/Section2All_xls.xlsx', 'T', '-A', { '20100' } ) ];
SourceData = [ SourceData; ObtainSingleBEASpreadsheet( 'https://apps.bea.gov/national/Release/XLS/Survey/Section3All_xls.xlsx', 'T', '-A', { '30100' } ) ];
SourceData = [ SourceData; ObtainSingleBEASpreadsheet( 'https://apps.bea.gov/national/Release/XLS/Survey/Section6All_xls.xlsx', 'T', '-A', { '60900B', '60900C', '60900D' } ) ];

save SourceData SourceData;
