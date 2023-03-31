function [newtracks, splits3] = newmergelinker(events, tracks,...
    splits, points, areas, intensities, t_cutoff)
%events{k} is an array where is row corresponds to either a merge or split
%and is given by
% events{k}(i,:) = [merge_track, other track, time, merge/split]
merge = 0;
split = 1;
idxCol = 1;
tCol = 2;
%Go through and for each event in each central branch, break up the track
%into subtracks
num_merges = sum(cellfun(@(x) sum(x(:,4) == merge), events));
merge_links = zeros(2*num_merges, 4);

nevents = cellfun(@(x) size(x,1), events);
ntime = size(tracks{1},1);
subtracks = cell(size(tracks,1), max(nevents) + 1);
trackidxs = zeros(length(tracks),1);
nstracks = ones(length(tracks),1);
for k = 1:length(events)
    nstracks(events{k}(1,1)) = nevents(k) + 1; %number of subtracks is equal to number of events+1
    for j = 1:nevents(k)
        %break each track with n splits or merges into n+1 subtracks
        tridx = events{k}(j,1); %track index
        etime = events{k}(j,3); %event time
        if j == 1
            prevtime = 0;
        else
            prevtime = events{k}(j-1,3);
        end
        times = transpose(prevtime+1:etime);
        %each subtrack has a column for indices of the point in each frame
        %and a column for time points
        subtrack = [tracks{tridx}(times), times];
%         if isempty(subtrack(~isnan(subtrack(:,1)),:))
%             disp(strcat("Subtrack (", num2str(tridx), ",", num2str(j), ") is empty"))
%             types = ["merge","split"];
%             disp(strcat("Previous event type is ", types(events{k}(j-1,end)+1)))
%             disp(strcat("Next event type is ", types(events{k}(j,end)+1)))
%         end
        subtracks{tridx,j} = subtrack(~isnan(subtrack(:,1)),:); %discard NaN rows
        
        if j == nevents(k) %add one more subtrack after the last event
            last_times = transpose((etime+1):ntime);
            subtrack = [tracks{tridx}(last_times), last_times];
            subtracks{tridx,j+1} = subtrack(~isnan(subtrack(:,1)),:);
        end
        trackidxs(tridx) = 1;
    end
end

%make subtracks for all other tracks (not in events)
otracks = find(trackidxs == 0);
for ti = 1:length(otracks)
    tridx = otracks(ti);
    subtrack = [tracks{tridx}, (1:ntime)'];
    subtracks{tridx,1} = subtrack(~isnan(subtrack(:,1)),:);
end

%  nstracks = sum(cellfun(@isempty,subtracks),2); %number of subtracks for each track (column vector)
% nstracks = sum(~cellfun(@isempty,subtracks),2); %number of subtracks for each track (column vector)

%keep track of which tracks are tacked to the end of existing tracks 
%instead of being their own track
kept_tracks = ones(size(subtracks,1),1);

count = 1;
%iterate through events and assign subtracks at merges and splits
disp('Linking merged subtracks')
for k = 1:length(events)
    if mod(ti, round(0.05*length(events))) == 0
        fprintf('.')
    end
    for si = 1:nevents(k)
        %use flag to differentiate between when we are linking to splitting
        %tracks (true) and when we are linking to the mid track (false)
        flag = false;
        eventType = events{k}(si,4); %split or merge (1 or 0)?
        if eventType == merge
            tridx = events{k}(si,1); %track index
            otidx = events{k}(si,2); %other track index
            etime = events{k}(si,3); %event time
            
            %get time, index, x, y, intensity, and area for input 1
            in1 = [tridx, si];
            if isempty(subtracks{in1(1),in1(2)})
                in1 = [tridx, si-1];
            end
            ii1 = subtracks{in1(1),in1(2)}(end, idxCol);
            ti1 = subtracks{in1(1),in1(2)}(end, tCol);
            xyi1 = points{ti1}(ii1,:);
            Ii1 = intensities{ti1}(ii1,:);
            Ai1 = areas{ti1}(ii1,:);
            
            %get time, index, x, y, intensity, and area for input 2
            in2 = [otidx, nstracks(otidx)];
            ii2 = subtracks{in2(1),in2(2)}(end, idxCol);
            ti2 = subtracks{in2(1),in2(2)}(end, tCol);
            xyi2 = points{ti2}(ii2,:);
            Ii2 = intensities{ti2}(ii2,:);
            Ai2 = areas{ti2}(ii2,:);
            
            %determine whether to link to mid or split track 
            if si < nevents(k) %if there is another event after this one
                neventType = events{k}(si+1,4); %next event type
                netime = events{k}(si+1,3); %next event time
                %get cell area just before splitting
                mid = [tridx, si+1];
                if ~isempty(subtracks{mid(1), mid(2)})
                    midx = subtracks{mid(1), mid(2)}(end, idxCol);
                    mtime = subtracks{mid(1), mid(2)}(end, tCol);
                    AbS = areas{mtime}(midx);
                else
                    %this deals with instances in which the same point on a
                    %given track had both a merge to it and a split from it
                    AbS = Ai1 + Ai2;
                end
                % if next event is a split within a given cutoff of the
                % previous merge and does not appear to be real cell
                % division (use relative size of nucleus before merge and just before split)
                if (neventType == split) && (netime <= etime + t_cutoff) && (AbS > max(Ai1, Ai2))
                    flag = true;
                    out1 = [events{k}(si+1,1), si+2]; %next segment of this track
                    io1 = subtracks{out1(1),out1(2)}(1, idxCol);
                    to1 = subtracks{out1(1),out1(2)}(1, tCol);
                    out2 = [events{k}(si+1,2), 1]; %start of other track
                    io2 = subtracks{out2(1),out2(2)}(1, idxCol);
                    to2 = subtracks{out2(1),out2(2)}(1, tCol);
                end
            end
            
            %Do the linking
            if flag %assume cells overlap and then diverge
                %get x, y, area, intensity for both output subtracks
                xyo1 = points{to1}(io1,:);
                Io1 = intensities{to1}(io1);
                Ao1 = areas{to1}(io1);
                xyo2 = points{to2}(io2,:);
                Io2 = intensities{to2}(io2);
                Ao2 = areas{to2}(io2);
                %make relative costs based on intensity, area, position
                %link in1 to out 1
                c11 = sum((xyi1 - xyo1).^2)*(1 + 2*abs(Ii1-Io1)/(Ii1+Io1))*(1 + 2*abs(Ai1-Ao1)/(Ai1+Ao1));
                %link in1 to out2
                c12 = sum((xyi1 - xyo2).^2)*(1 + 2*abs(Ii1-Io2)/(Ii1+Io2))*(1 + 2*abs(Ai1-Ao2)/(Ai1+Ao2));
                %link in2 to out1
                c21 = sum((xyi2 - xyo1).^2)*(1 + 2*abs(Ii2-Io1)/(Ii2+Io1))*(1 + 2*abs(Ai2-Ao1)/(Ai2+Ao1));
                %link in2 to out 2
                c22 = sum((xyi2 - xyo2).^2)*(1 + 2*abs(Ii2-Io2)/(Ii2+Io2))*(1 + 2*abs(Ai2-Ao2)/(Ai2+Ao2));
                %if c11 + c22 <= c12 + c21, link in1->out1, in2->out2
                %otherwise linke in1->out2, in2->out1
                if c11 + c22 <= c12 + c21
                    merge_links(count, :) = [in1, out1];
                    merge_links(count + 1, :) = [in2, out2];
                else
                    merge_links(count, :) = [in1, out2];
                    merge_links(count + 1, :) = [in2, out1];
                end
                count = count + 2;
                %remove this split from the list of splitting events - it
                %is now essentially a double gap-close
                discard_split = [events{k}(si+1,2), events{k}(si+1,1)];
                splits = splits(~ismember(splits(:,1:2), discard_split, 'rows'), :);
%                 kept_tracks(out1(1)) = 0;
                kept_tracks(out2(2)) = 0;
            else %assume one cell disappeared or died
                %get time, index, x, y, intensity, and area for start of mid track
                out1 = [tridx, si+1];
                io1 = subtracks{out1(1),out1(2)}(1, idxCol);
                to1 = subtracks{out1(1),out1(2)}(1, tCol);
                xyo1 = points{to1}(io1,:);
                Io1 = intensities{to1}(io1,:);
                Ao1 = areas{to1}(io1,:);
                %make costs and stuff
                c1 = sum((xyi1 - xyo1).^2)*(1 + 2*abs(Ii1-Io1)/(Ii1+Io1))*(1 + 2*abs(Ai1-Ao1)/(Ai1+Ao1));
                c2 = sum((xyi2 - xyo1).^2)*(1 + 2*abs(Ii2-Io1)/(Ii2+Io1))*(1 + 2*abs(Ai2-Ao1)/(Ai2+Ao1));
                if c1 >= c2
                    merge_links(count,:) = [in1, out1];
                    merge_links(count+1,:) = [in2, 0, 0];
                else
                    merge_links(count,:) = [in2, out1];
                    merge_links(count+1,:) = [in1, 0, 0];
                end
                count = count + 2;
            end
        end
    end
end
fprintf('\n')

% %%%%%% may need to add a mid_track field column to merge_links to cut out 
% discontinuous segments of tracks %%%%%%

%go through subtracks and make lists of continuously strung together
%subtracks
%merge_links(k,:) = [start_track, subtrack idx, end_track, subtrack idx]
disp('Reorganizing subtracks within tracks')
subtracklists = cell(size(subtracks, 1),1);
scount = 1;
for ti = 1:size(subtracks,1)
    if mod(ti, round(0.05*size(subtracks,1))) == 0
        fprintf('.')
    end
    if kept_tracks(ti)
        %current track
        ctrack = ti;
        %current subtrack idx
        cstrack = 1;
        %number of subtracks in current track
        ctracklen = nstracks(ctrack);
        %initialize the list of subtracks
        stracklist = [];
        criterion = true;
        %starting at first subtrack of current track, keep adding subtracks
        %sequentially; add all subtracks up until this track is an input for
        %another track; then switch the current track and subtrack to that
        %output (without failing to add intermediary subtrack if necessary -
        %if there was a merge then split assignment, both inputs get the
        %subtrack in between the merge and the split)
        loop_count = 0;
        while criterion
            %find links for which this track is an input track at a node
            ins = merge_links(ismember(merge_links(:,1),ctrack),:);
            %if there are not any merges with this track as an input, exit the
            %loop at the end of this iteration
            if isempty(ins)
                criterion =  false;
                last_strack = ctracklen;
            else
                %only consider links at or after the current subtrack
                ins = ins(ins(:,2) >= cstrack,:);
                if isempty(ins)
                    %if no more merges from this track, add the remaining
                    %subtrack idxs for this track to the list and exit the loop
                    %at the end of the iteration
                    last_strack = ctracklen;
                    criterion = false;
                else
                    [last_strack, lidx] = min(ins(:,2)); %subtrack idx of next merge from this track
                    next_strack = ins(lidx,3:4); %track and subtrack idxs of output of next merge
                end

            end
            %add subtracks from current to last subtrack before next merge to
            %the list
            stracklist = [stracklist;ctrack*ones(last_strack-cstrack+1,1),(cstrack:last_strack)'];
            if criterion
                ctrack = next_strack(1);
                cstrack = next_strack(2);
                if ctrack == 0
                    criterion = false;
                else
                    ctracklen = nstracks(ctrack);
                end
            end
            loop_count = loop_count + 1;
            if loop_count > 50
                error('infinite loop')
            end
        end
        subtracklists{scount} = stracklist;
        scount = scount + 1;
    end
end
fprintf('\n')

subtracklists = subtracklists(1:(scount-1)); %discard empty entries at the bottom

%rearrange the old subtracks using the new assignments
%adjust the track indices in splits accordingly for splits3
newtracks = cell(size(subtracklists,1),1);
for ti = 1:size(subtracklists,1)
    newtracks{ti} = NaN(ntime,1);
    for si = 1:size(subtracklists{ti},1)
        ctrack = subtracklists{ti}(si,1);
        cstrack = subtracklists{ti}(si,2);
        subtrack = subtracks{ctrack, cstrack};
        newtracks{ti}(subtrack(:,2)) = subtrack(:,1);
    end
end

%adjust the indices of splits for new assignments
%splits(i,:) = [start_track, mid_track, split_time, cost]
%we want to iterate through splits and replace the old track indices with
%new indices based on where the subtracks were shuffled to
splits3 = splits;
for si = 1:size(splits,1)
    midx = splits(si, 2); %midtrack idx
    sidx = splits(si, 1); %start track idx
    stime = splits(si, 3); %time of split (frame in time-lapse)
    %subtrack for stidx is first subtrack of track stidx:
    s_st = [sidx, 1];
    %find appropriate subtrack for split in mid track
    for ti = 1:nstracks(midx)
        if ismember(stime, subtracks{midx, ti}(:,tCol))
            m_st = [midx, ti];
        end
    end
    %find which tracks each of these subtracks is now in
    %this probably gives unfavorable behavior if duplicate subtrack
    %assignments are not properly sorted out
    for ti = 1:size(subtracklists,1)
        %replace track index of split track in splits
        if ismember(s_st, subtracklists{ti}, 'rows')
            splits3(si,1) = ti;
        end
        %replace index of mid track in splits
        if ismember(m_st, subtracklists{ti}, 'rows')
            splits3(si,2) = ti;
        end
    end
end

% error('need to examine some of these variables')







end