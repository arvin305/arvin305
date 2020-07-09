function [TradeRestrictionBoolean,IND_flipped,IND] = IndicatorDataPreconditioning(TradeRestrictionBoolean,IND,AC)
%% This function preconditions the data for an indicator
    TradeRestrictionBooleanFlipped = flipud(TradeRestrictionBoolean);
    %IND positions to Carry Forward
    IND0 = IND;
    IND_flipped = flipud(IND);
    IND_CFpos = isnan(IND) - isnan(AC);
    IND_CFpos_flippedUD = flipud(IND_CFpos);
    for i = 1:size(IND,2)
        if size(find(IND_CFpos_flippedUD(:,i) ~= 1,1, 'first'),1) > 0
            IND_FirstNonNaNPos(i) = find(IND_CFpos_flippedUD(:,i) ~= 1,1, 'first');
            if IND_FirstNonNaNPos(i) ~= 1
                NaNPos = find(IND_CFpos_flippedUD(:,i) == 1);
                NaNPos = NaNPos(NaNPos > IND_FirstNonNaNPos(i));
                for k = 1:length(NaNPos)
                    IND_flipped(NaNPos(k),i) = IND_flipped(NaNPos(k)-1,i);
                end
                TradeRestrictionBooleanFlipped(1:IND_FirstNonNaNPos(i)-1,i) = TradeRestrictionBooleanFlipped(1:IND_FirstNonNaNPos(i)-1,i) + 1;
            else if sum(find(IND_CFpos_flippedUD(:,i) ~= 1)) ~= sum((1:size(IND,1)))
                    NaNPos = find(IND_CFpos_flippedUD(:,i) == 1);
                    for k = 1:length(NaNPos)
                        IND_flipped(NaNPos(k),i) = IND_flipped(NaNPos(k)-1,i);
                    end
                end
            end
        else
            IND_FirstNonNaNPos(i) = 0;
            TradeRestrictionBooleanFlipped(:,i) = TradeRestrictionBooleanFlipped(:,i) + 1;
        end
% % Uncomment lines below to follow Initial and Post-conditioning Columns next to Adjusted Close (AC) Indicator
%         IND(:,i) = flipud(IND_flipped(:,i));
%         [IND0(:,i),IND(:,i),AC(:,i)]   
    end
    TradeRestrictionBooleanFlipped(isnan(IND_flipped)) = 1;
    IND = flipud(IND_flipped);
    TradeRestrictionBoolean = flipud(TradeRestrictionBooleanFlipped);
    TradeRestrictionBoolean(TradeRestrictionBoolean > 1) = 1;
end

% [IND_flipped,isnan(IND_flipped),TradeRestrictionBooleanFlipped]