function Table = PrefaceVariableNames( Table, Preface )

    Table.Properties.VariableNames = cellfun( @( x ) [ Preface x ], Table.Properties.VariableNames, 'UniformOutput', false );

end
