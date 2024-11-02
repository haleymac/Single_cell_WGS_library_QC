#Dictionary Creation Functions
import os 


def get_output_extensions():
    #Get an dictionary of output file extensions.
    output_extensions = {'output': {'wgs': '.wgs_metrics.txt',
                                     'gc_bias': '.gc_bias_metrics.txt',
                                     'gc_bias_summary': '.gc_bias_summary.txt',
                                     'gc_bias_plot': '.gc_bias_metrics_chart.pdf',
                                     'insert': '.insert_metrics.txt',
                                     'insert_plot': '.insert_metrics_histo.pdf',
                                     'duplicate': '.duplicate_metrics.txt',
                                     'mt_coverage': '.mt_coverage.bed',
                                     'rna' : '.RNA_Metrics.tsv'}}

    return output_extensions




def create_output_file_dict(cell_ids, output_dir):
    #Create output dictionaries for each cell and file type.
    extensions = get_output_extensions()

    file_dict = {}

    for file_category in extensions.keys():
        file_dict[file_category] = {}

        for file_type, file_extension in extensions[file_category].items():
            file_dict[file_category][file_type] = {}

            for cell_id in cell_ids:
                file_dict[file_category][file_type][cell_id] = os.path.join(output_dir, "individual_cell_metrics",
                                                                            '{}{}'.format(cell_id, file_extension))

    return file_dict


def create_input_bam_dict(cell_ids, input_dir):
    #Create a dictionary of single cell BAM files given a subset of cell IDs.
    bam_files = [x for x in os.listdir(input_dir) if x.endswith('.bam')]

    bam_paths = [os.path.join(input_dir, bam_file) for bam_file in bam_files]

    bam_cell_ids = [x.split('.')[0].split('_')[0] for x in bam_files]

    full_bam_dict = dict(zip(bam_cell_ids, bam_paths))

    bam_dict = dict((cell_id, full_bam_dict[cell_id]) for cell_id in cell_ids)

    return bam_dict



def create_input_bam_dict(cell_ids, input_dir):
    #Create a dictionary of single cell BAM files given a subset of cell IDs.
    bam_files = [x for x in os.listdir(input_dir) if x.endswith('.bam')]

    bam_paths = [os.path.join(input_dir, bam_file) for bam_file in bam_files]

    bam_cell_ids = [x.split('.')[0] for x in bam_files]

    full_bam_dict = dict(zip(bam_cell_ids, bam_paths))

    bam_dict = dict((cell_id, full_bam_dict[cell_id]) for cell_id in cell_ids)

    return bam_dict