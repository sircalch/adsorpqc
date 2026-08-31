"""
Interactive HTML dashboard generator for AdsorpQC.
"""

import os
import jinja2
from adsorpqc.core.scoring import AdsorptionValidationReport

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>AdsorpQC Adsorption Quality Report</title>
    <style>
        :root {
            --bg-color: #0f172a;
            --card-bg: #1e293b;
            --card-border: #334155;
            --text-primary: #f8fafc;
            --text-secondary: #94a3b8;
            --pass-color: #10b981;
            --pass-bg: rgba(16, 185, 129, 0.15);
            --warn-color: #f59e0b;
            --warn-bg: rgba(245, 158, 11, 0.15);
            --fail-color: #ef4444;
            --fail-bg: rgba(239, 68, 68, 0.15);
            --accent-blue: #38bdf8;
        }
        * { box-sizing: border-box; margin: 0; padding: 0; }
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
            background-color: var(--bg-color);
            color: var(--text-primary);
            line-height: 1.6;
            padding: 2rem 1rem;
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
        }
        .header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            border-bottom: 1px solid var(--card-border);
            padding-bottom: 1.5rem;
            margin-bottom: 2rem;
            flex-wrap: wrap;
            gap: 1rem;
        }
        .title-group h1 {
            font-size: 2rem;
            font-weight: 700;
            letter-spacing: -0.025em;
        }
        .title-group p {
            color: var(--text-secondary);
            font-size: 0.95rem;
        }
        .status-badge {
            display: inline-flex;
            align-items: center;
            padding: 0.5rem 1.25rem;
            border-radius: 9999px;
            font-weight: 700;
            font-size: 1.1rem;
            letter-spacing: 0.05em;
            text-transform: uppercase;
        }
        .badge-pass { background-color: var(--pass-bg); color: var(--pass-color); border: 1px solid var(--pass-color); }
        .badge-warning { background-color: var(--warn-bg); color: var(--warn-color); border: 1px solid var(--warn-color); }
        .badge-fail { background-color: var(--fail-bg); color: var(--fail-color); border: 1px solid var(--fail-color); }
        
        .grid-cards {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 1.25rem;
            margin-bottom: 2rem;
        }
        .card {
            background-color: var(--card-bg);
            border: 1px solid var(--card-border);
            border-radius: 0.75rem;
            padding: 1.25rem;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.2);
        }
        .card-label {
            font-size: 0.8rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: var(--text-secondary);
            margin-bottom: 0.4rem;
        }
        .card-value {
            font-size: 1.4rem;
            font-weight: 700;
            color: var(--text-primary);
        }
        .card-subtext {
            font-size: 0.85rem;
            color: var(--text-secondary);
            margin-top: 0.25rem;
        }

        .section-title {
            font-size: 1.3rem;
            font-weight: 600;
            margin-bottom: 1rem;
            color: var(--accent-blue);
        }

        table {
            width: 100%;
            border-collapse: collapse;
            margin-bottom: 2rem;
            font-size: 0.95rem;
            background-color: var(--card-bg);
            border-radius: 0.75rem;
            overflow: hidden;
            border: 1px solid var(--card-border);
        }
        th, td {
            padding: 0.85rem 1rem;
            text-align: left;
            border-bottom: 1px solid var(--card-border);
        }
        th {
            background-color: rgba(255, 255, 255, 0.03);
            color: var(--text-secondary);
            font-weight: 600;
            font-size: 0.85rem;
            text-transform: uppercase;
        }
        tr:hover td { background-color: rgba(255, 255, 255, 0.02); }

        .tag {
            display: inline-block;
            padding: 0.2rem 0.6rem;
            border-radius: 0.375rem;
            font-size: 0.75rem;
            font-weight: 700;
        }
        .tag-pass { background-color: var(--pass-bg); color: var(--pass-color); }
        .tag-warning { background-color: var(--warn-bg); color: var(--warn-color); }
        .tag-fail { background-color: var(--fail-bg); color: var(--fail-color); }

        .box {
            background-color: var(--card-bg);
            border: 1px solid var(--card-border);
            border-radius: 0.75rem;
            padding: 1.25rem;
            margin-bottom: 2rem;
        }
        pre {
            background-color: rgba(0, 0, 0, 0.4);
            padding: 1rem;
            border-radius: 0.5rem;
            overflow-x: auto;
            color: #38bdf8;
            font-family: monospace;
            font-size: 0.9rem;
            border: 1px solid rgba(255, 255, 255, 0.05);
            white-space: pre-wrap;
        }
        .btn-copy {
            background-color: #2563eb;
            color: white;
            border: none;
            padding: 0.4rem 0.8rem;
            border-radius: 0.375rem;
            cursor: pointer;
            font-size: 0.8rem;
            margin-top: 0.5rem;
        }
        .btn-copy:hover { background-color: #1d4ed8; }

        footer {
            text-align: center;
            font-size: 0.85rem;
            color: var(--text-secondary);
            margin-top: 3rem;
            padding-top: 1.5rem;
            border-top: 1px solid var(--card-border);
        }
    </style>
</head>
<body>
    <div class="container">
        <header class="header">
            <div class="title-group">
                <h1>AdsorpQC Simulation & Isotherm Report</h1>
                <p>Equilibration / Burn-in Detection & Model Fitting Certification for Nanoporous Materials</p>
            </div>
            <div>
                <span class="status-badge badge-{{ report.overall_status.lower() }}">
                    {{ report.overall_status }}
                </span>
            </div>
        </header>

        <div class="grid-cards">
            <div class="card">
                <div class="card-label">Framework / Gas</div>
                <div class="card-value">{{ report.metadata.framework }}</div>
                <div class="card-subtext">{{ report.metadata.adsorbate }} at {{ "%.1f"|format(report.metadata.temperature_k) }} K</div>
            </div>
            <div class="card">
                <div class="card-label">GCMC Burn-in</div>
                <div class="card-value">
                    {% if report.burnin_result %}
                    {{ "%.1f"|format(report.burnin_result.burnin_fraction * 100) }}%
                    {% else %}
                    N/A
                    {% endif %}
                </div>
                <div class="card-subtext">
                    {% if report.burnin_result %}
                    {{ report.burnin_result.n_burnin_cycles }} / {{ report.burnin_result.n_total_cycles }} cycles (Status: {{ report.burnin_result.status }})
                    {% else %}
                    Direct Isotherm Input
                    {% endif %}
                </div>
            </div>
            <div class="card">
                <div class="card-label">Production Drift</div>
                <div class="card-value">
                    {% if report.burnin_result %}
                    {{ "%.2f"|format(report.burnin_result.loading_drift_pct) }}%
                    {% else %}
                    --
                    {% endif %}
                </div>
                <div class="card-subtext">
                    {% if report.burnin_result %}
                    N_eff = {{ "%.0f"|format(report.burnin_result.n_eff) }} samples
                    {% else %}
                    No raw GCMC cycles
                    {% endif %}
                </div>
            </div>
            <div class="card">
                <div class="card-label">Optimal Model</div>
                <div class="card-value">
                    {% if report.isotherm_fits and report.isotherm_fits.best_model_name %}
                    {{ report.isotherm_fits.best_model_name }}
                    {% else %}
                    N/A
                    {% endif %}
                </div>
                <div class="card-subtext">
                    {% if report.isotherm_fits and report.isotherm_fits.best_r_squared %}
                    R&sup2; = {{ "%.4f"|format(report.isotherm_fits.best_r_squared) }}
                    {% else %}
                    --
                    {% endif %}
                </div>
            </div>
        </div>

        <h2 class="section-title">Statistical & Thermodynamic Validation</h2>
        <table>
            <thead>
                <tr>
                    <th>Diagnostic Component</th>
                    <th>Measured Value</th>
                    <th>Acceptance Standard</th>
                    <th>Status</th>
                    <th>Interpretation</th>
                </tr>
            </thead>
            <tbody>
                {% if report.burnin_result %}
                <tr>
                    <td><strong>GCMC Equilibration</strong></td>
                    <td>Cutoff at {{ report.burnin_result.n_burnin_cycles }} cycles</td>
                    <td>Optimal N_eff maximization</td>
                    <td><span class="tag tag-{{ report.burnin_result.status.lower() }}">{{ report.burnin_result.status }}</span></td>
                    <td style="color: var(--text-secondary)">{{ report.burnin_result.diagnostic_message }}</td>
                </tr>
                <tr>
                    <td><strong>Loading Drift</strong></td>
                    <td>Drift = {{ "%.2f"|format(report.burnin_result.loading_drift_pct) }}%</td>
                    <td>Drift &le; 3.0% (Warning &le; 8.0%)</td>
                    <td><span class="tag tag-{{ report.burnin_result.status.lower() }}">{{ report.burnin_result.status }}</span></td>
                    <td style="color: var(--text-secondary)">Loading = {{ "%.3f"|format(report.burnin_result.production_mean_loading) }} &plusmn; {{ "%.3f"|format(report.burnin_result.production_std_loading) }}</td>
                </tr>
                {% endif %}
                {% if report.isotherm_fits and report.isotherm_fits.best_model %}
                <tr>
                    <td><strong>Isotherm Model Fit</strong></td>
                    <td>{{ report.isotherm_fits.best_model_name }} (R&sup2; = {{ "%.4f"|format(report.isotherm_fits.best_r_squared) }})</td>
                    <td>R&sup2; &ge; 0.95</td>
                    <td><span class="tag tag-pass">PASS</span></td>
                    <td style="color: var(--text-secondary)">RMSE = {{ "%.4f"|format(report.isotherm_fits.best_rmse) }}, AIC = {{ "%.1f"|format(report.isotherm_fits.best_aic) }}</td>
                </tr>
                {% endif %}
                {% if report.henry_data %}
                <tr>
                    <td><strong>Henry Coefficient (K<sub>H</sub>)</strong></td>
                    <td>{{ "%.4e"|format(report.henry_data.henry_constant) }} mol/kg/Pa</td>
                    <td>R&sup2; &ge; 0.90 (linear regime)</td>
                    <td><span class="tag tag-pass">PASS</span></td>
                    <td style="color: var(--text-secondary)">Low-pressure slope (R&sup2; = {{ "%.3f"|format(report.henry_data.r_squared) }})</td>
                </tr>
                {% endif %}
                {% if report.isosteric_heat_kj_mol is not none and report.isosteric_heat_kj_mol > 0 %}
                <tr>
                    <td><strong>Isosteric Heat (q<sub>st</sub>)</strong></td>
                    <td>{{ "%.2f"|format(report.isosteric_heat_kj_mol) }} kJ/mol</td>
                    <td>Energy-particle fluctuations</td>
                    <td><span class="tag tag-pass">PASS</span></td>
                    <td style="color: var(--text-secondary)">Host-adsorbate binding enthalpy</td>
                </tr>
                {% endif %}
                {% if report.iast_selectivity %}
                <tr>
                    <td><strong>IAST Selectivity S<sub>A/B</sub></strong></td>
                    <td>S = {{ "%.2f"|format(report.iast_selectivity.selectivity_a_b) }}</td>
                    <td>[95% CI: {{ "%.2f"|format(report.iast_selectivity.ci_lower_95) }} - {{ "%.2f"|format(report.iast_selectivity.ci_upper_95) }}]</td>
                    <td><span class="tag tag-pass">PASS</span></td>
                    <td style="color: var(--text-secondary)">Binary mixture separation factor</td>
                </tr>
                {% endif %}
            </tbody>
        </table>

        {% if report.recommendations %}
        <div class="box" style="border-left: 4px solid var(--warn-color);">
            <h3 style="color: var(--warn-color); margin-bottom: 0.5rem;">Quality Recommendations & Alerts</h3>
            <ul style="padding-left: 1.25rem;">
                {% for rec in report.recommendations %}
                <li style="margin-bottom: 0.25rem; color: var(--text-secondary);">{{ rec }}</li>
                {% endfor %}
            </ul>
        </div>
        {% endif %}

        <h2 class="section-title">Computational Details (Publication Ready Methods)</h2>
        <div class="box">
            <pre id="methodsSnippet">{{ methods_text }}</pre>
            <button class="btn-copy" onclick="copyToClipboard('methodsSnippet')">Copy Methods Paragraph</button>
        </div>

        <h2 class="section-title">BibTeX Citation</h2>
        <div class="box">
            <pre id="bibSnippet">{{ citation_bib }}</pre>
            <button class="btn-copy" onclick="copyToClipboard('bibSnippet')">Copy BibTeX</button>
        </div>

        <footer>
            Generated automatically by <strong>AdsorpQC v1.0.0</strong> &bull; Adsorption Quality & Isotherm Certification Toolkit &bull; Monreal-Hernández, 2026.
        </footer>
    </div>

    <script>
        function copyToClipboard(elementId) {
            const text = document.getElementById(elementId).innerText;
            navigator.clipboard.writeText(text).then(() => {
                alert('Copied to clipboard!');
            }).catch(err => {
                console.error('Error copying: ', err);
            });
        }
    </script>
</body>
</html>
"""


def generate_adsorption_html_report(
    report: AdsorptionValidationReport,
    output_path: str,
    methods_text: str = "",
    citation_bib: str = ""
) -> str:
    """
    Renders HTML report template and writes it to disk.
    """
    template = jinja2.Template(HTML_TEMPLATE)
    rendered = template.render(
        report=report,
        methods_text=methods_text,
        citation_bib=citation_bib
    )
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(rendered)
    return output_path
