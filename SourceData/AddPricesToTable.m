function Table = AddPricesToTable( Table, Variables )

    for i = 1 : length( Variables )
        
        Variable = Variables{ i };
        
        Table.( [ 'Price' Variable ] ) = Table.( [ 'Nominal' Variable ] ) ./ Table.( [ 'Real' Variable ] );
        
    end
    
end
