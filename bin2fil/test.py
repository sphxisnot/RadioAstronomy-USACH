from main import file_runthrough, write_header
from sampling_qol import load_data
from convertidor import convert_and_write

test_file= load_data("tstFold04.params")
write_header(test_file)
convert_and_write(test_file)
