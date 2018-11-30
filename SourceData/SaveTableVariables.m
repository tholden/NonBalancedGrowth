function SaveTableVariables( Table, FileName )

    Struct = table2struct( Table, 'ToScalar', true ); %#ok<NASGU>
    
    save( FileName, '-struct', 'Struct' );

end
