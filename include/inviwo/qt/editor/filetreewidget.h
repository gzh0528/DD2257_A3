/*********************************************************************************
 *
 * Inviwo - Interactive Visualization Workshop
 *
 * Copyright (c) 2019 Inviwo Foundation
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *********************************************************************************/
#pragma once

#include <inviwo/qt/editor/inviwoqteditordefine.h>
#include <inviwo/core/common/inviwo.h>

#include <warn/push>
#include <warn/ignore/all>

#include <QTreeView>
#include <QAbstractItemModel>
#include <QIcon>

#include <warn/pop>

class QSortFilterProxyModel;

namespace inviwo {

class InviwoApplication;
class TreeItem;
class TreeModel;

class IVW_QTEDITOR_API FileTreeWidget : public QTreeView {
#include <warn/push>
#include <warn/ignore/all>
    Q_OBJECT
#include <warn/pop>
public:
    enum ListElemType { File = 1, Section, SubSection, None };

    enum ItemRoles { FileName = Qt::UserRole + 100, Path, Type, ExampleWorkspace };

    explicit FileTreeWidget(InviwoApplication* app, QWidget* parent = nullptr);
    virtual ~FileTreeWidget() = default;

    void updateRecentWorkspaces(const QStringList& recentFiles);
    void updateExampleEntries();
    void updateRegressionTestEntries();

    bool selectRecentWorkspace(int index);

    void setFilter(const QString& str);

signals:
    void selectedFileChanged(const QString& filename, bool isExample);
    void loadFile(const QString& filename, bool isExample);

private:
    InviwoApplication* inviwoApp_;

    TreeModel* model_;
    QSortFilterProxyModel *proxyModel_;

    TreeItem* recentWorkspaceItem_ = nullptr;
    TreeItem* examplesItem_ = nullptr;
    TreeItem* regressionTestsItem_ = nullptr;

    QIcon fileIcon_;
};

class IVW_QTEDITOR_API TreeItem {
public:
    explicit TreeItem(TreeItem* parent = nullptr);
    TreeItem(const QString& caption, FileTreeWidget::ListElemType type, TreeItem* parent = nullptr);
    TreeItem(const QIcon& icon, const std::string& filename, bool isExample = false,
             TreeItem* parent = nullptr);
    virtual ~TreeItem();

    void addChild(TreeItem* child);
    void addChildren(std::vector<TreeItem*> children);

    bool insertChildren(int position, int count);
    bool removeChildren(int position, int count);
    void removeChildren();

    TreeItem* child(int row);
    int row() const;
    int childCount() const;
    int columnCount() const;
    TreeItem* parent() const;

    virtual QVariant data(int column, int role) const;
    virtual FileTreeWidget::ListElemType type() const;

    void setData(const QString& caption, FileTreeWidget::ListElemType type);
    void setData(const QIcon& icon, const std::string& filename, bool isExample);

private:
    TreeItem* parent_;
    std::vector<TreeItem*> childItems_;

    FileTreeWidget::ListElemType type_ = FileTreeWidget::ListElemType::None;

    QIcon icon_;
    QString caption_;
    QString file_;
    QString path_;
    bool isExample_;
};

class IVW_QTEDITOR_API TreeModel : public QAbstractItemModel {
#include <warn/push>
#include <warn/ignore/all>
    Q_OBJECT
#include <warn/pop>
public:
    explicit TreeModel(QObject* parent = nullptr);
    ~TreeModel();

    virtual QModelIndex index(int row, int column,
                              const QModelIndex& parent = QModelIndex()) const override;
    virtual Qt::ItemFlags flags(const QModelIndex& index) const override;
    virtual QVariant data(const QModelIndex& index, int role) const override;
    virtual QVariant headerData(int section, Qt::Orientation orientation,
                                int role = Qt::DisplayRole) const override;
    virtual QModelIndex parent(const QModelIndex& index) const override;
    virtual int rowCount(const QModelIndex& parent = QModelIndex()) const override;
    virtual int columnCount(const QModelIndex& parent = QModelIndex()) const override;

    virtual bool insertRows(int position, int rows,
                            const QModelIndex& parent = QModelIndex()) override;
    virtual bool removeRows(int position, int rows,
                            const QModelIndex& parent = QModelIndex()) override;

    void updateCategory(TreeItem* item, std::vector<TreeItem*> children);

    void addEntry(TreeItem* root, TreeItem* child);
    bool removeEntry(TreeItem* node);
    bool removeChildren(TreeItem* root);
    QModelIndex getIndex(TreeItem* item, int column = 0) const;

private:
    TreeItem* getItem(const QModelIndex& index) const;

    TreeItem* root_;
};

}  // namespace inviwo
