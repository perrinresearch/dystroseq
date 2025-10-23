"""
Generate reproducibility reports for exonSeq analysis runs.
"""

import json
import os
import sys
import platform
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional

from . import VERSION_INFO


def get_system_info() -> Dict[str, str]:
    """Get system information for reproducibility."""
    return {
        "platform": platform.platform(),
        "python_version": sys.version,
        "python_executable": sys.executable,
        "architecture": platform.architecture()[0],
        "machine": platform.machine(),
        "processor": platform.processor(),
        "hostname": platform.node()
    }


def get_git_info() -> Dict[str, str]:
    """Get git information if available."""
    try:
        # Get git commit hash
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"], 
            capture_output=True, 
            text=True, 
            cwd=Path(__file__).parent.parent
        )
        commit_hash = result.stdout.strip() if result.returncode == 0 else "unknown"
        
        # Get git branch
        result = subprocess.run(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"], 
            capture_output=True, 
            text=True, 
            cwd=Path(__file__).parent.parent
        )
        branch = result.stdout.strip() if result.returncode == 0 else "unknown"
        
        return {
            "commit_hash": commit_hash,
            "branch": branch,
            "repository": "exonSeq"
        }
    except (subprocess.SubprocessError, FileNotFoundError):
        return {
            "commit_hash": "unknown",
            "branch": "unknown", 
            "repository": "exonSeq"
        }


def get_dependencies_info() -> Dict[str, str]:
    """Get information about installed dependencies."""
    try:
        # Try to get pip list output
        result = subprocess.run(
            [sys.executable, "-m", "pip", "list", "--format=json"], 
            capture_output=True, 
            text=True
        )
        if result.returncode == 0:
            packages = json.loads(result.stdout)
            # Filter for relevant packages
            relevant_packages = {}
            for pkg in packages:
                name = pkg.get("name", "").lower()
                if any(keyword in name for keyword in ["requests", "biopython", "numpy", "pandas"]):
                    relevant_packages[pkg["name"]] = pkg["version"]
            return relevant_packages
    except (subprocess.SubprocessError, json.JSONDecodeError):
        pass
    
    return {"note": "Could not determine package versions"}


def generate_reproducibility_report(
    command_args: Dict[str, Any],
    output_dir: Path,
    analysis_type: str = "export"
) -> Path:
    """
    Generate a comprehensive reproducibility report.
    
    Args:
        command_args: Dictionary of command line arguments used
        output_dir: Output directory where report should be saved
        analysis_type: Type of analysis performed (export, exon, etc.)
    
    Returns:
        Path to the generated report file
    """
    timestamp = datetime.now().isoformat()
    
    # Gather all information
    report_data = {
        "exonSeq_info": VERSION_INFO,
        "timestamp": timestamp,
        "analysis_type": analysis_type,
        "command_arguments": command_args,
        "system_info": get_system_info(),
        "git_info": get_git_info(),
        "dependencies": get_dependencies_info(),
        "output_directory": str(output_dir.absolute()),
        "reproducibility": {
            "note": "This report contains all information needed to reproduce this analysis",
            "exonSeq_version": VERSION_INFO["version"],
            "command_to_reproduce": _generate_reproduce_command(command_args, analysis_type)
        }
    }
    
    # Create report filename
    report_filename = f"exonSeq_report_{analysis_type}_{timestamp.replace(':', '-').split('.')[0]}.json"
    report_path = output_dir / report_filename
    
    # Write JSON report
    with open(report_path, 'w', encoding='utf-8') as f:
        json.dump(report_data, f, indent=2, ensure_ascii=False)
    
    # Also create a human-readable summary
    summary_path = output_dir / f"exonSeq_summary_{analysis_type}_{timestamp.replace(':', '-').split('.')[0]}.txt"
    _write_human_readable_summary(report_data, summary_path)
    
    return report_path


def _generate_reproduce_command(command_args: Dict[str, Any], analysis_type: str) -> str:
    """Generate a command string that can be used to reproduce the analysis."""
    if analysis_type == "export":
        cmd_parts = ["python", "-m", "exonseq", "export"]
        
        # Add required arguments
        if "species" in command_args:
            cmd_parts.extend(["--species", str(command_args["species"])])
        if "gene" in command_args:
            cmd_parts.extend(["--gene", str(command_args["gene"])])
        if "reference" in command_args:
            cmd_parts.extend(["--reference", str(command_args["reference"])])
        if "coords" in command_args:
            cmd_parts.extend(["--coords", str(command_args["coords"])])
        if "transcript" in command_args:
            cmd_parts.extend(["--transcript", str(command_args["transcript"])])
        if "variant" in command_args:
            cmd_parts.extend(["--variant", str(command_args["variant"])])
        if "out" in command_args:
            cmd_parts.extend(["--out", str(command_args["out"])])
            
    elif analysis_type == "exon":
        cmd_parts = ["python", "-m", "exonseq", "exon"]
        # Add exon-specific arguments
        if "exon_command" in command_args:
            cmd_parts.append(command_args["exon_command"])
        if "exon_args" in command_args:
            cmd_parts.extend(command_args["exon_args"])
    else:
        cmd_parts = ["python", "-m", "exonseq", analysis_type]
    
    return " ".join(cmd_parts)


def _write_human_readable_summary(report_data: Dict[str, Any], summary_path: Path) -> None:
    """Write a human-readable summary of the analysis."""
    with open(summary_path, 'w', encoding='utf-8') as f:
        f.write("=" * 60 + "\n")
        f.write("exonSeq ANALYSIS REPORT\n")
        f.write("=" * 60 + "\n\n")
        
        f.write(f"Analysis Type: {report_data['analysis_type']}\n")
        f.write(f"Timestamp: {report_data['timestamp']}\n")
        f.write(f"exonSeq Version: {report_data['exonSeq_info']['version']}\n\n")
        
        f.write("COMMAND TO REPRODUCE:\n")
        f.write("-" * 30 + "\n")
        f.write(f"{report_data['reproducibility']['command_to_reproduce']}\n\n")
        
        f.write("COMMAND ARGUMENTS:\n")
        f.write("-" * 30 + "\n")
        for key, value in report_data['command_arguments'].items():
            f.write(f"{key}: {value}\n")
        f.write("\n")
        
        f.write("SYSTEM INFORMATION:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Platform: {report_data['system_info']['platform']}\n")
        f.write(f"Python Version: {report_data['system_info']['python_version']}\n")
        f.write(f"Architecture: {report_data['system_info']['architecture']}\n")
        f.write(f"Hostname: {report_data['system_info']['hostname']}\n\n")
        
        f.write("GIT INFORMATION:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Repository: {report_data['git_info']['repository']}\n")
        f.write(f"Branch: {report_data['git_info']['branch']}\n")
        f.write(f"Commit: {report_data['git_info']['commit_hash']}\n\n")
        
        f.write("OUTPUT DIRECTORY:\n")
        f.write("-" * 30 + "\n")
        f.write(f"{report_data['output_directory']}\n\n")
        
        f.write("REPRODUCIBILITY NOTE:\n")
        f.write("-" * 30 + "\n")
        f.write("This analysis was performed with exonSeq. To reproduce these results,\n")
        f.write("use the command shown above with the same version of exonSeq and\n")
        f.write("compatible system environment.\n")


def get_version_string() -> str:
    """Get a formatted version string for display."""
    return f"exonSeq v{VERSION_INFO['version']}"


def print_version_info() -> None:
    """Print version information to stdout."""
    print(f"exonSeq version {VERSION_INFO['version']}")
    print(f"Author: {VERSION_INFO['author']}")
    print(f"Description: {VERSION_INFO['description']}")
