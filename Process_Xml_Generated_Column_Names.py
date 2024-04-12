import sys
from collections import Counter

def refine_segments(fields):
    # Lowercase and split fields
    fields = [field.lower() for field in fields]
    split_fields = [field.split('_') for field in fields]

    # Flatten the split list to count occurrences of each segment
    all_segments = [segment for sublist in split_fields for segment in sublist]
    segment_count = Counter(all_segments)
    total_lists = len(fields)

    # Identify segments that appear only once and those that are common
    unique_segments = {segment for segment, count in segment_count.items() if count == 1}
    common_segments = {segment for segment, count in segment_count.items() if count >= 0.1 * total_lists}

    # Set to track all modified fields to ensure uniqueness
    seen_fields = set()
    output_fields = []

    # Process the fields to remove common segments
    for split_field in split_fields:
        filtered_segments = [segment for segment in split_field if segment not in common_segments]

        if not filtered_segments and split_field:
            # If removing common segments would empty the list, use the least common segment
            least_common_segment = min(split_field, key=lambda x: segment_count[x])
            filtered_segments.append(least_common_segment)
        
        new_field = '_'.join(filtered_segments)

        if new_field not in seen_fields:
            seen_fields.add(new_field)
            output_fields.append(new_field)
        else:
            # If duplicate, add the least common unique segment to make it unique
            unique_addition = min((seg for seg in split_field if seg not in seen_fields), 
                                  key=lambda x: segment_count[x], default='')
            if unique_addition:
                filtered_segments.append(unique_addition)
                new_field = '_'.join(filtered_segments)
                seen_fields.add(new_field)
                output_fields.append(new_field)

    return output_fields

def main():
    # Read from stdin, assume each input line is a single string of tab-delimited fields
    input_lines = [line.strip() for line in sys.stdin if line.strip()]
    fields = [item for line in input_lines for item in line.split('\t')]

    # Process the fields to refine and make them unique
    refined_fields = refine_segments(fields)

    # Output the final unique fields
    for field in refined_fields:
        print(field)

if __name__ == "__main__":
    main()
