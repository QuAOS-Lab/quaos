import argparse
import os


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--r', type=str, required=True, help='Path to requirements'
    )
    return parser.parse_args()


def find_requirements_file(path):
    if os.path.isfile(path):
        return path
    elif os.path.isfile(os.path.join(os.getcwd(), path)):
        return os.path.join(os.getcwd(), path)
    else:
        raise FileNotFoundError(f"Requirements file not found: {path}")


if __name__ == '__main__':
    args = parse_args()
    requirements_file = find_requirements_file(args.r)
    requirements_file = os.path.abspath(requirements_file).replace('\\', '/')
    output_dir = os.path.dirname(requirements_file)

    # Read current requirements
    with open(output_dir + '/current_requirements.txt', 'r') as f:
        current_requirements = {line.split('==')[0]: line.strip() for line in f}

    # Read existing requirements
    with open(requirements_file, 'r') as f:
        existing_requirements = {
            line.split('==')[0]: line.strip() for line in f
        }

    # Update requirements
    updated_requirements = []
    for package, version in existing_requirements.items():
        if package in current_requirements:
            updated_requirements.append(current_requirements[package])
        else:
            updated_requirements.append(version)

    # Write updated requirements
    with open(output_dir + '/updated_requirements.txt', 'w') as f:
        for requirement in updated_requirements:
            f.writelines(requirement + '\n')
