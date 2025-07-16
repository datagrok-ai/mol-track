from data_generator import DataGenerator


if __name__ == "__main__":
    dg = DataGenerator()

    dg.generate_schema_and_mapping_files(output_dir="./tests/performance_tests/data")
