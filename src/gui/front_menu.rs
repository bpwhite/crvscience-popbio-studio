// Copyright (C) 2026 CRV Science Network (https://crvscience.com)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://gnu.org>.

use egui::{Color32, Margin, RichText, Ui};

use crate::model::{ModuleGroup, ModuleId};
use crate::workset::{SaveSlotInfo, WorksetState};

pub enum FrontMenuAction {
    None,
    LaunchModule(ModuleId),
    SaveNamedSlot(String),
    LoadNamedSlot(String),
    DeleteNamedSlot(String),
    ExportWorksetFile,
    ImportWorksetFile,
}

/// Renders the organized Front Menu ("Module Hub") grouping all 19 population biology,
/// theoretical genetics, and empirical sequence/haplotype modules into 5 taxonomic domains,
/// alongside the Workset ("Save Game") slot manager.
pub fn render_front_menu(
    ui: &mut Ui,
    workset: &mut WorksetState,
    filter_query: &mut String,
    new_slot_name: &mut String,
    slots: &[SaveSlotInfo],
) -> FrontMenuAction {
    let mut action = FrontMenuAction::None;

    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            egui::Frame::group(ui.style())
                .fill(Color32::from_rgb(18, 26, 38))
                .inner_margin(Margin::same(14))
                .show(ui, |ui| {
                    ui.horizontal_wrapped(|ui| {
                        ui.vertical(|ui| {
                            ui.label(
                                RichText::new(
                                    "PopBio Studio — Population Biology, Evolutionary Genetics & Haplotype Network Workstation",
                                )
                                .size(20.0)
                                .strong()
                                .color(Color32::from_rgb(236, 244, 255)),
                            );
                            ui.label(
                                RichText::new(
                                    "Select any organismal, theoretical, or empirical sequence module below, or resume a saved Workset checkpoint.",
                                )
                                .color(Color32::from_rgb(168, 184, 204)),
                            );
                        });
                    });

                    ui.add_space(8.0);
                    ui.horizontal(|ui| {
                        ui.label(RichText::new("Workset Title:").strong());
                        ui.text_edit_singleline(&mut workset.workset_title);
                        ui.separator();
                        ui.label("Filter Modules:");
                        ui.text_edit_singleline(filter_query);
                        if !filter_query.is_empty() && ui.button("Clear").clicked() {
                            filter_query.clear();
                        }
                    });
                });

            ui.add_space(10.0);

            egui::CollapsingHeader::new(
                RichText::new(
                    "Workset Checkpoints & 'Save Game' Slot Manager (SQLite + .popworkset.json)",
                )
                .size(15.0)
                .strong(),
            )
            .default_open(true)
            .show(ui, |ui| {
                ui.horizontal_wrapped(|ui| {
                    ui.label("Slot Name:");
                    ui.text_edit_singleline(new_slot_name);
                    if ui.button("Save Checkpoint Slot").clicked() {
                        let clean = if new_slot_name.trim().is_empty() {
                            "Slot_1".to_string()
                        } else {
                            new_slot_name.trim().to_string()
                        };
                        action = FrontMenuAction::SaveNamedSlot(clean);
                    }
                    if ui.button("Quicksave (F5)").clicked() {
                        action = FrontMenuAction::SaveNamedSlot("Quicksave".to_string());
                    }
                    if ui.button("Export .popworkset.json...").clicked() {
                        action = FrontMenuAction::ExportWorksetFile;
                    }
                    if ui.button("Load .popworkset.json...").clicked() {
                        action = FrontMenuAction::ImportWorksetFile;
                    }
                });

                ui.add_space(4.0);
                ui.horizontal(|ui| {
                    ui.label("Study Notes:");
                    ui.text_edit_singleline(&mut workset.study_notes);
                });

                if !slots.is_empty() {
                    ui.add_space(6.0);
                    egui::Grid::new("saved_workset_slots_grid")
                        .striped(true)
                        .min_col_width(110.0)
                        .show(ui, |ui| {
                            ui.label(RichText::new("Slot").strong());
                            ui.label(RichText::new("Workset Title").strong());
                            ui.label(RichText::new("Active Module").strong());
                            ui.label(RichText::new("Saved (UTC)").strong());
                            ui.label(RichText::new("Actions").strong());
                            ui.end_row();

                            for slot in slots {
                                ui.label(RichText::new(&slot.slot_name).monospace());
                                ui.label(&slot.workset_title);
                                ui.label(&slot.active_module);
                                ui.label(slot.saved_at_utc.chars().take(19).collect::<String>());
                                ui.horizontal(|ui| {
                                    if ui.button("Resume").clicked() {
                                        action =
                                            FrontMenuAction::LoadNamedSlot(slot.slot_name.clone());
                                    }
                                    if ui.small_button("Delete").clicked() {
                                        action = FrontMenuAction::DeleteNamedSlot(
                                            slot.slot_name.clone(),
                                        );
                                    }
                                });
                                ui.end_row();
                            }
                        });
                } else {
                    ui.label(
                        RichText::new(
                            "No saved SQLite slots yet — press 'Quicksave (F5)' or 'Save Checkpoint Slot' anytime to snapshot all 19 modules.",
                        )
                        .italics()
                        .color(Color32::GRAY),
                    );
                }
            });

            ui.add_space(12.0);

            let query = filter_query.trim().to_ascii_lowercase();
            for group in ModuleGroup::ALL {
                let [r, g, b] = group.accent_rgb();
                let accent = Color32::from_rgb(r, g, b);

                let group_modules: Vec<ModuleId> = ModuleId::ALL
                    .into_iter()
                    .filter(|m| m.group() == group)
                    .filter(|m| {
                        if query.is_empty() {
                            true
                        } else {
                            let meta = m.metadata();
                            meta.title.to_ascii_lowercase().contains(&query)
                                || meta.organisms.to_ascii_lowercase().contains(&query)
                                || meta.summary.to_ascii_lowercase().contains(&query)
                                || meta.code.to_ascii_lowercase().contains(&query)
                        }
                    })
                    .collect();

                if group_modules.is_empty() {
                    continue;
                }

                egui::Frame::group(ui.style())
                    .fill(Color32::from_rgb(20, 24, 32))
                    .inner_margin(Margin::same(12))
                    .show(ui, |ui| {
                        ui.horizontal(|ui| {
                            ui.label(
                                RichText::new(format!(" {} ", group.short_badge()))
                                    .strong()
                                    .color(Color32::WHITE)
                                    .background_color(accent),
                            );
                            ui.label(
                                RichText::new(group.title())
                                    .size(16.5)
                                    .strong()
                                    .color(Color32::from_rgb(235, 240, 248)),
                            );
                        });
                        ui.label(
                            RichText::new(group.subtitle())
                                .size(12.5)
                                .color(Color32::from_rgb(160, 174, 192)),
                        );
                        ui.add_space(8.0);

                        ui.columns(group_modules.len().clamp(1, 3), |cols| {
                            for (idx, &mod_id) in group_modules.iter().enumerate() {
                                let col_ui = &mut cols[idx % cols.len()];
                                let meta = mod_id.metadata();
                                let phase = workset.current_phase(mod_id);
                                let [pr, pg, pb] = phase.badge_rgb();

                                egui::Frame::group(col_ui.style())
                                    .fill(Color32::from_rgb(28, 34, 46))
                                    .inner_margin(Margin::same(10))
                                    .show(col_ui, |card| {
                                        card.horizontal(|card| {
                                            card.label(
                                                RichText::new(format!("[{}]", meta.code))
                                                    .monospace()
                                                    .strong()
                                                    .color(accent),
                                            );
                                            card.label(
                                                RichText::new(meta.title)
                                                    .strong()
                                                    .color(Color32::WHITE),
                                            );
                                        });
                                        card.label(
                                            RichText::new(meta.organisms)
                                                .italics()
                                                .size(12.0)
                                                .color(Color32::from_rgb(138, 201, 255)),
                                        );
                                        card.add_space(4.0);
                                        card.label(
                                            RichText::new(meta.summary)
                                                .size(12.0)
                                                .color(Color32::from_rgb(200, 208, 220)),
                                        );
                                        card.add_space(6.0);
                                        card.horizontal_wrapped(|card| {
                                            card.label(
                                                RichText::new(format!(" {} ", phase.label()))
                                                    .size(11.0)
                                                    .color(Color32::WHITE)
                                                    .background_color(Color32::from_rgb(pr, pg, pb)),
                                            );
                                            if card.button("Open Module ->").clicked() {
                                                action = FrontMenuAction::LaunchModule(mod_id);
                                            }
                                        });
                                    });
                                col_ui.add_space(6.0);
                            }
                        });
                    });
                ui.add_space(10.0);
            }
        });

    action
}
