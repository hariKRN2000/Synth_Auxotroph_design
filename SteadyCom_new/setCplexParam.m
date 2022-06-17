function LP = setCplexParam(LP, solverParam,verbFlag)
%Set the parameters of the CPLEX object according to the structure solverParam
%For example:
%[solverParam.simplex.display, solverParam.tune.display, solverParam.barrier.display,...
%     solverParam.sifting.display, solverParam.conflict.display] = deal(0);
%[solverParam.simplex.tolerances.optimality, solverParam.simplex.tolerances.feasibility] = deal(1e-9,1e-8);
if nargin < 3
    verbFlag = true;
end
if isempty(fieldnames(solverParam))
    return
end
[paramList, paramPath] = getParamList(LP.Param, 0);
[paramUserList, paramUserPath] = getParamList(solverParam, 1);
paramIden = false(numel(paramUserList), 1);
for p = 1:numel(paramUserList)
    f = strcmpi(paramList,paramUserList{p});
    if sum(f) == 1
        paramIden(p) = true;
        str = ['LP.Param.' paramPath{f} '.Cur = solverParam.' paramUserPath{p} ';'];
        eval(str);
    elseif sum(f) > 1
        if ismember(lower(paramUserPath{p}), paramPath);
            paramIden(p) = true;
            str = ['LP.Param.' lower(paramUserPath{p}) '.Cur = solverParam.' paramUserPath{p} ';'];
            eval(str);
        else
            if verbFlag
                fprintf('solverParam.%s cannot be uniquely identified as a valid cplex parameter. Ignore.\n', paramUserPath{p});
            end
        end
    else
        if verbFlag
            fprintf('solverParam.%s cannot be identified as a valid cplex parameter. Ignore.\n', paramUserPath{p});
        end
    end
end
end

function [paramList, paramPath] = getParamList(param, bottomFlag)
%for matching CPLEX parameters appropriately
structCur = param;
lv = 1;
lvFieldN = zeros(10,1);
lvFieldN(1) = 1;
lvField = cell(10, 1);
lvField{lv} = fieldnames(structCur);
paramPath = {};
paramList = {};
while lv > 0
    if isstruct(structCur.(lvField{lv}{lvFieldN(lv)})) 
        if ~isempty(fieldnames(structCur.(lvField{lv}{lvFieldN(lv)})))
            structCur = structCur.(lvField{lv}{lvFieldN(lv)});
            lv = lv + 1;
            lvFieldN(lv) = 1;
            lvField{lv} = fieldnames(structCur);
        else
            while lvFieldN(lv) == numel(lvField{lv})
                lv = lv - 1;
                if lv == 0
                    break
                end
            end
            if lv > 0
                lvFieldN(lv) = lvFieldN(lv) + 1;
                structCur = param;
                for j = 1:lv-1
                    structCur = structCur.(lvField{j}{lvFieldN(j)});
                end
            end
        end
    else
        if ~bottomFlag
            lv = lv - 1;
        end
        if lv > 0
            c = {};
            for j = 1:lv
                c = [c lvField{j}(lvFieldN(j))];
            end
            paramPath = [paramPath; strjoin(c,'.')];
            paramList = [paramList; c(end)];
            while lvFieldN(lv) == numel(lvField{lv})
                lv = lv - 1;
                if lv == 0
                    break
                end
            end
            if lv > 0
                lvFieldN(lv) = lvFieldN(lv) + 1;
                structCur = param;
                for j = 1:lv-1
                    structCur = structCur.(lvField{j}{lvFieldN(j)});
                end
            end
        else
            lv = 1;
            if lvFieldN(lv) == numel(lvField{lv})
                break
            else
                lvFieldN(1) = lvFieldN(1) + 1;
            end
        end
    end
end

end