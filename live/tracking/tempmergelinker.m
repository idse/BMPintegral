function [merge_links, discard_splits] = tempmergelinker(events, t_cutoff)
%Output of this function should indicate whether each merged track
%terminates or links to a split track, and, if it splits, which track it
%links to; this is equivalent to an additional gap closing step
%Each row should be as follows:
%merge_links = [merged_track, split_track] ~> [end_track, start_track]
%make this [merge_track, split_track, mid_track, t_merge, t_split]
%If split_track = 0, then merged_track terminates; otherwise, it is
%gap-closed with split_track

%Add an output to show which splits have been converted to gap-closes
%discard_splits(i,:) = [split_track, mid_track]

%find the total number of merges
num_merges = sum(cellfun(@(x) sum(x(:,4) == 0), events));
num_splits = sum(cellfun(@(x) sum(x(:,4) == 1), events));
merge_links = zeros(num_merges, 5);
discard_splits = zeros(num_splits, 2);
%use count to index the merge_links array
count = 1;
count2 = 1;

for k = 1:length(events)
    for j = 1:(size(events{k},1)-1)
        
        if events{k}(j,4) == 0
            in1 = events{k}(j,1);
            in2 = events{k}(j,2);
            time = events{k}(j,3);
            time_next = events{k}(j+1,3);
            next_event = events{k}(j+1,4);
            if (next_event == 1) && (time_next - time <= t_cutoff)
                out2 = events{k}(j+1, 2);
                merge_links(count,:) = [in2, out2, in1, time, time_next];
                discard_splits(count2,:) = [out2, in1];
                count2 = count2 + 1;
            else
                merge_links(count,:) = [in2, zeros(1,4)];
            end
            count = count + 1;
        end
    end
    if events{k}(end,4) == 0
        in2 = events{k}(end,2);
        merge_links(count,:) = [in2, zeros(1,4)];
        count = count + 1;
    end
end

discard_splits = discard_splits(discard_splits(:,1) ~= 0,:);

end