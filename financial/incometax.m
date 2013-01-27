function incomeOUT = incometax(Income)
 % Tax     Rates    Band    Diff
 FedTax  = [0.15,   41544, 41544
            0.22,   83088, 50712
            0.26,   133800, 0];
        
 ProvTax = [0.0506,  36146, 36147
            0.077,   72293, 10708
            0.105,   83001, 17786
            0.1229,  100787, 0];

for p = 1:2
    
    % Fed Tax Deductions
    Index  = Income(p) > FedTax(:,2);
    if sum(Index) == 0
        Deduct = Income(p)*FedTax(1,1);
    else
    over   = find(Index,1,'last');
    Deduct = sum(FedTax(Index,3).*FedTax(Index,1)) +...
        (Income(p) - FedTax(over,2))*FedTax(over+1,1);
    end
    % Prov Tax Deductions
    Index  = Income(p) > ProvTax(:,2);
    if sum(Index) == 0
        Deduct = Deduct + Income(p)*ProvTax(1,1);
    else
    over   = find(Index,1,'last');
    Deduct = Deduct + sum(ProvTax(Index,3).*ProvTax(Index,1)) +...
        (Income(p) - ProvTax(over,2))*ProvTax(over+1,1);
    end
    % Income after Deductions
    incomeOUT(p) = Income(p) - Deduct;
end
end

