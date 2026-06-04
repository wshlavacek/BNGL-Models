from pathlib import Path

import yaml


def test_repo_skills_have_required_frontmatter() -> None:
    skills_dir = Path("skills")
    skill_files = sorted(skills_dir.rglob("SKILL.md"))

    assert skill_files

    for skill_file in skill_files:
        content = skill_file.read_text(encoding="utf-8")
        assert content.startswith("---\n"), f"{skill_file} missing YAML frontmatter"

        _, frontmatter_text, _ = content.split("---", 2)
        frontmatter = yaml.safe_load(frontmatter_text)

        assert isinstance(frontmatter, dict), f"{skill_file} frontmatter is not a mapping"
        assert frontmatter.get("name"), f"{skill_file} missing frontmatter name"
        assert frontmatter.get("description"), f"{skill_file} missing frontmatter description"
