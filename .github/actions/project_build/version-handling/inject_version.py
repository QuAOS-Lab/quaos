#!/usr/bin/env python3
import re
import sys
import argparse


def generate_version(fallback_version="0.0.0.0a0", include_distance=True):
    try:
        import setuptools_scm

        def custom_version_scheme(version):
            return version.format_with("{tag}a{distance}")

        def custom_local_scheme(version):
            return ""

        version = setuptools_scm.get_version(
            version_scheme=custom_version_scheme,
            local_scheme=custom_local_scheme,
            fallback_version=fallback_version,
            tag_regex=r'^v?(\d+\.\d+\.\d+)(?:\.\d+)?$'
        )

        # If include_distance is False, strip the 'aN' suffix
        if not include_distance and 'a' in version:
            version = version.split('a')[0]

        return version
    except Exception:
        return fallback_version


def inject_version(version, pyproject_path="pyproject.toml"):
    try:
        with open(pyproject_path, 'r', encoding='utf-8') as f:
            content = f.read()
    except FileNotFoundError:
        return False

    content = re.sub(
        r'dynamic\s*=\s*\["version"\]',
        f'version = "{version}"',
        content
    )

    try:
        with open(pyproject_path, 'w', encoding='utf-8') as f:
            f.write(content)
        return True
    except Exception:
        return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fallback-version', default='0.0.0.0a0')
    parser.add_argument('--pyproject-path', default='pyproject.toml')
    parser.add_argument('--include-distance', default='true', choices=['true', 'false'],
                        help='Include distance (aN) in version for non-exact tags')

    args = parser.parse_args()

    include_distance = args.include_distance.lower() == 'true'

    version = generate_version(
        fallback_version=args.fallback_version,
        include_distance=include_distance
    )

    success = inject_version(version, args.pyproject_path)

    print(f"version={version}")
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
